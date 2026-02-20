% AfdmSystem: AFDM通信链路仿真系统类
classdef AfdmSystem < handle

    properties
        % --- 基础参数 ---
        BitsPerSymbol
        ModulationOrder
        NumDataSubcarriers          % 有效数据的子载波数
        NumPaths
        
        % --- 帧结构参数 ---
        MaxNormDoppler
        MaxPathDelays
        PrefixLength
        TotalSubcarriers
        
        % --- AFDM 核心参数 ---
        DopplerGuard
        ChirpParam1
        ChirpParam2
        ZeroPaddingLength
        
        % --- 导频与有效数据 --- 
        PilotIndex
        PilotSymbol
        PilotSnr          %  导频信噪比
        NumActiveCarriers
        ActiveIndices

        % --- 波形与均衡器配置 ---
        WaveformType            % 波形类型
        CsiMode                 % CSI模式 (Perfect, Estimated)
        EqualizerType           % 均衡器类型
        NumMaxIterations        % DFE最大迭代次数

    end

    methods(Access = public)
        % --- 构造函数 ---
        function obj = AfdmSystem()
            % 这里设置默认参数，也可以改为输入参数
            obj.BitsPerSymbol = 2;
            obj.ModulationOrder = 2^obj.BitsPerSymbol;
            obj.NumDataSubcarriers = 2^8;
            obj.NumPaths = 3;
            obj.MaxNormDoppler = 2;
            % obj.PathDelays = [2 3 5];
            obj.MaxPathDelays = 5;
            obj.DopplerGuard = 4;       % k_v
            obj.PilotIndex = 1;
            obj.PilotSymbol = 1 + 1j;
            obj.PilotSnr = 35;          % 导频信噪比 (dB)
            
            updateDerivedParams(obj);
        end

        % --- 计算衍生参数 ---
        function updateDerivedParams(obj)
            % maxDelay = max(obj.PathDelays);
            obj.PrefixLength = obj.MaxPathDelays;
            obj.TotalSubcarriers = obj.NumDataSubcarriers + obj.PrefixLength;
            
            % 计算 c1, c2
            obj.ChirpParam1 = (2 * (obj.MaxNormDoppler + obj.DopplerGuard) + 1) / (2 * obj.NumDataSubcarriers);
            obj.ChirpParam2 = 1 / (obj.NumDataSubcarriers^2 * 2 * pi);
            
            if (2 * (obj.MaxNormDoppler + obj.DopplerGuard) * (obj.MaxPathDelays + 1)) + obj.MaxPathDelays > obj.NumDataSubcarriers      %必须满足正交条件
                error("子载波不满足正交！\n");
            end

            % 计算 ZP 和 有效索引
            obj.ZeroPaddingLength = (obj.MaxPathDelays + 1) * (2 * (obj.MaxNormDoppler + obj.DopplerGuard) + 1) - 1;
            
            % 计算有效载荷
            obj.NumActiveCarriers = obj.NumDataSubcarriers - 2 * obj.ZeroPaddingLength - 1;
            obj.ActiveIndices = (obj.ZeroPaddingLength + 2) : (obj.ZeroPaddingLength + 1 + obj.NumActiveCarriers);
            
            % 简单的检查
            if max(obj.ActiveIndices) > (obj.NumDataSubcarriers - obj.ZeroPaddingLength)
                error('参数配置错误: ZP太长或子载波太少');
            end
        end

        % --- 发送机 ---
        function [txSignal, originalData] = transmit(obj)
            % 生成数据
            originalData = randi([0 obj.ModulationOrder-1], obj.NumActiveCarriers, 1);
            qamSymbols = qammod(originalData, obj.ModulationOrder, 'UnitAveragePower', true);
            
            % 组帧
            % 全0初始化
            oneFrame = zeros(obj.NumDataSubcarriers, 1);

            if obj.CsiMode == "Estimated"
                oneFrame(obj.PilotIndex) = obj.PilotSymbol; % 插入导频
            end

            oneFrame(obj.ActiveIndices) = qamSymbols;
            
            %  变换回时域 IDAFT / IDFT
            if upper(obj.WaveformType) == "AFDM"
                timeFrame = runInverseDaftTransform(obj, oneFrame);
            else
                timeFrame = runInverseDftTransform(obj,oneFrame);
            end
            
            %  前缀 CP / CPP
            if upper(obj.WaveformType) == "AFDM" 
                    preFix = timeFrame(end - obj.PrefixLength + 1 : end) .* ...
                    exp(-1j * 2 * pi * obj.ChirpParam1 * (obj.NumDataSubcarriers^2 + 2 * obj.NumDataSubcarriers * (-obj.PrefixLength:-1).'));
            else
                    preFix = timeFrame(end - obj.PrefixLength + 1 : end);
            end
            txSignal = [preFix; timeFrame];
        end

        % --- 接收机 ---
        function rxData = receive(obj, rxSignal, noisePower,physicalChannelMatrix)
            % 移除Prefix
            rxNoPrefix = rxSignal(obj.PrefixLength + 1 : end);
            
            % DAFT / DFT
            if upper(obj.WaveformType) == "AFDM"
                finalRx = runDaftTransform(obj, rxNoPrefix);
            else
                finalRx = runDftTransform(obj,rxNoPrefix);
            end
            % daftRx = runDaftTransform(obj, rxNoPrefix);
            
            if obj.CsiMode == "Estimated"
                % 提取导频进行信道估计
                [effectiveChannelMatrix, ~] = runChannelEstimator(obj, finalRx);
            else
                % 上帝视角 (对比AFDM OFDM)
                effectiveChannelMatrix = generateEffectiveChannelMatrix(obj,physicalChannelMatrix);
            end

            
            % 均衡
            eqSignal = equalizer(obj,finalRx, effectiveChannelMatrix, noisePower);
            
            % 解调
            rxData = qamdemod(eqSignal, obj.ModulationOrder, 'UnitAveragePower', true);
            rxData = rxData(:);                                 % 确保是列向量
        end

        % --- 信道估计与重建 ---
        % 仅AFDM系统可用
        function [estimatedEffectiveChannelMatrix, finalEstimatedParams] = runChannelEstimator(obj, receivedSignalDaft)
            % --- 配置搜索范围 ---
            delaySearchRange = 0 : obj.PrefixLength;              
            dopplerSearchRange = -obj.MaxNormDoppler : obj.MaxNormDoppler; 
            
            % --- 初始化 ---
            residualSignal = receivedSignalDaft; % 残差信号
            estimatedPaths = zeros(obj.NumPaths, 3); % [Delay, Doppler, Gain]
            
            % 窗半径 (根据你的 ZP 长度调整，通常 2-4 足够)
            % 目的：只保留导频主峰附近的能量，切断远处拖尾
            winRadius = 3; 

            for p = 1 : obj.NumPaths
                % 1. 粗搜索 (Coarse Search)
                maxCorr = -inf;
                coarseDelay = 0;
                coarseIntDoppler = 0;
                
                for dIdx = 1 : length(delaySearchRange)
                    currDelay = delaySearchRange(dIdx);
                    for dopIdx = 1 : length(dopplerSearchRange)
                        currDop = dopplerSearchRange(dopIdx);
                        
                        % 生成探测向量
                        probeVector = buildDaftResponse(obj, currDelay, currDop);
                        
                        % 计算相关性
                        correlation = abs(probeVector' * residualSignal);
                        
                        if correlation > maxCorr
                            maxCorr = correlation;
                            coarseDelay = currDelay;
                            coarseIntDoppler = currDop;
                        end
                    end
                end
                
                % 2. 精搜索 (Fine Search)
                costFunc = @(fracDop) - abs(buildDaftResponse(obj, coarseDelay, coarseIntDoppler + fracDop)' * residualSignal);
                
                [bestFracDoppler, ~] = fminbnd(costFunc, -0.5, 0.5, optimset('TolX', 1e-6));
                finalDoppler = coarseIntDoppler + bestFracDoppler;
                
                % 3. 生成全长基向量 (未加窗)
                fullBasisVector = buildDaftResponse(obj, coarseDelay, finalDoppler);

                % 4. === 关键步骤：构造动态窗函数 ===
                % 自动寻找基向量的峰值位置作为窗中心 (因为 AFDM 的峰值位置取决于 c1 和 delay，计算复杂，直接找峰值最准)
                [~, peakIndex] = max(abs(fullBasisVector)); 
                
                % 生成循环索引 (处理边界情况)
                windowIndicesRaw = (peakIndex - winRadius : peakIndex + winRadius).';
                % 转换为有效的 1-based 索引 (1 ~ N)
                windowIndices = mod(windowIndicesRaw - 1, obj.NumDataSubcarriers) + 1;
                
                % 创建掩膜 (Mask)
                mask = zeros(obj.NumDataSubcarriers, 1);
                mask(windowIndices) = 1; % 矩形窗
                
                % 应用窗函数：得到“纯净”的导频响应
                cleanBasisVector = fullBasisVector .* mask;

                % 5. === 基于纯净响应重新估计增益 ===
                % 使用最小二乘法，但只在窗内进行。这样远处的噪声不会污染增益估计。
                % 公式: h = (p_clean' * r) / (p_clean' * p_clean) / pilot
                energy = cleanBasisVector' * cleanBasisVector;
                
                if energy > 1e-10 % 防止除零
                    estimatedGain = (cleanBasisVector' * residualSignal) / energy / obj.PilotSymbol;
                else
                    estimatedGain = 0;
                end
                
                % 6. 串行干扰消除 (SIC)
                % 注意：这里我们减去的是 [estimatedGain * cleanBasisVector]
                % 意味着我们只减去了主峰。远处的拖尾（可能是噪声，也可能是泄露）保留在残差中，
                % 避免了因模型不匹配导致的“错误减法”。
                pathComponent = estimatedGain * obj.PilotSymbol * cleanBasisVector;
                residualSignal = residualSignal - pathComponent;
                
                % 7. 存储结果
                estimatedPaths(p, :) = [coarseDelay, finalDoppler, estimatedGain];
            end
            
            % --- 全局最小二乘优化 (可选，保持原样) ---
            finalBasisMatrix = zeros(obj.NumDataSubcarriers, obj.NumPaths);
            for p = 1 : obj.NumPaths
                finalBasisMatrix(:, p) = buildDaftResponse(obj, estimatedPaths(p, 1), estimatedPaths(p, 2));
            end
            refinedGains = pinv(finalBasisMatrix * obj.PilotSymbol) * receivedSignalDaft;
            estimatedPaths(:, 3) = refinedGains;

            % --- 输出 ---
            finalEstimatedParams = estimatedPaths;
            tempPhysicalChannelMatrix = rebuildLtvChannel(obj, finalEstimatedParams);
            estimatedEffectiveChannelMatrix = generateEffectiveChannelMatrix(obj, tempPhysicalChannelMatrix);
        end

        % --- 运行单次蒙特卡洛仿真 ---
        function ber = runMonteCarlo(obj, snrDb, targetErrors, maxFrames)
            noisePower = 10^(-snrDb / 10);
            pilotPower = 10^(obj.PilotSnr / 10) * noisePower;
            totalErrors = 0;
            totalBits = 0;
            frameCount = 0;
            
            fprintf('仿真开始 SNR: %d dB...\n', snrDb);
            
            while totalErrors < targetErrors && frameCount < maxFrames
                % 生成信道
                pathDelays = randperm(obj.MaxPathDelays + 1, obj.NumPaths) - 1;
                theta = (rand(1, obj.NumPaths) * 2 * pi) - pi;
                dopplers = obj.MaxNormDoppler * cos(theta);
                gains = (randn(1, obj.NumPaths) + 1j * randn(1, obj.NumPaths)) / sqrt(2 * obj.NumPaths);
                
                physicalChannelMatrix = LtvChannel(obj.TotalSubcarriers, pathDelays, dopplers, gains);
                % 发送前修正导频功率
                obj.PilotSymbol = sqrt(pilotPower) * (obj.PilotSymbol / abs(obj.PilotSymbol));
                % 发送
                [txSignal, txData] = transmit(obj);
                
                % 通过信道
                noise = sqrt(noisePower/2) * (randn(size(txSignal)) + 1j * randn(size(txSignal)));
                rxSignal = physicalChannelMatrix * txSignal + noise;
                
                % 接收
                rxData = receive(obj,rxSignal, noisePower,physicalChannelMatrix);
                
                % 统计
                [errs, ~] = biterr(txData, rxData);
                totalErrors = totalErrors + errs;
                totalBits = totalBits + obj.NumActiveCarriers * obj.BitsPerSymbol;
                frameCount = frameCount + 1;
            end
            
            ber = totalErrors / totalBits;
            fprintf('  -> 完成. 帧数: %d, BER: %.2e\n', frameCount, ber);
        end

    end

    methods (Access = private)
        % --- DAFT ---
        function demodulatedSignal = runDaftTransform(obj,receivedSignal)
            
            signalLength = size(receivedSignal, 1);

            chirpMatrix1 = diag(exp(-1j * 2 * pi * obj.ChirpParam1 * ((0:signalLength-1).^2)));
            chirpMatrix2 = diag(exp(-1j * 2 * pi * obj.ChirpParam2 * ((0:signalLength-1).^2)));
            dftMatrix = dftmtx(signalLength) ./ sqrt(signalLength);

            transformMatrix = chirpMatrix2 * dftMatrix * chirpMatrix1;
            demodulatedSignal = transformMatrix * receivedSignal;

        end

        % --- IDAFT ---
        function timeDomainSignal = runInverseDaftTransform(obj, daftDomainSignal)

            signalLength = size(daftDomainSignal, 1);

            chirpMatrix1 = diag(exp(-1j * 2 * pi * obj.ChirpParam1 * ((0:signalLength-1).^2)));
            chirpMatrix2 = diag(exp(-1j * 2 * pi * obj.ChirpParam2 * ((0:signalLength-1).^2)));
            dftMatrix = dftmtx(signalLength) ./ sqrt(signalLength);

            transformMatrix = chirpMatrix2 * dftMatrix * chirpMatrix1;
            % IDAFT 是 DAFT 的逆（对于酉矩阵是共轭转置）
            timeDomainSignal = transformMatrix' * daftDomainSignal;

        end

        % --- DFT ---
        function demodulatedSignal = runDftTransform(obj,receivedSignal)
            signalLength = size(receivedSignal, 1);
            dftMatrix = dftmtx(signalLength) ./ sqrt(signalLength);
            demodulatedSignal = dftMatrix * receivedSignal;
        end

        % --- IDFT ---
        function timeDomainSignal = runInverseDftTransform(obj,dftDomainSignal)
            signalLength = size(dftDomainSignal, 1);
            dftMatrix = dftmtx(signalLength) ./ sqrt(signalLength);
            timeDomainSignal = dftMatrix' * dftDomainSignal;
        end

        % --- 等效信道矩阵生成 ---
        function effectiveChannelMatrix = generateEffectiveChannelMatrix(obj,channelMatrixPhysical)
            % totalSubcarriers = size(channelMatrixPhysical, 1);
            % numDataSubcarriers = obj.TotalSubcarriers - obj.PrefixLength;

            % 这个索引包含了ZP部分
            dataIndices = (1 + obj.PrefixLength) : obj.TotalSubcarriers;
            
            % 构造 CPP 添加矩阵 (M)
            cppInsertionMatrix = zeros(obj.TotalSubcarriers, obj.NumDataSubcarriers);
            
            % 数据部分直接映射
            cppInsertionMatrix(obj.PrefixLength + 1 : end, :) = eye(obj.NumDataSubcarriers);
            
            if upper(obj.WaveformType) == "AFDM"
                % gamma_CPP
                gammaVector = exp(-1j * 2 * pi * obj.ChirpParam1 * (obj.NumDataSubcarriers^2 + 2 * obj.NumDataSubcarriers * (-obj.PrefixLength:-1).'));
            else
                gammaVector = ones(obj.PrefixLength, 1);
            end
            
            cppInsertionMatrix(1 : obj.PrefixLength, (obj.NumDataSubcarriers - obj.PrefixLength + 1) : obj.NumDataSubcarriers) = diag(gammaVector);
            % 计算等效时域矩阵 (N_data x N_data)
            % 截取接收信号的数据部分 = H(dataIndices, :) * M * S0
            effectiveChannelTime = channelMatrixPhysical(dataIndices, :) * cppInsertionMatrix;

            % 变换域矩阵
            dftMatrix  = dftmtx(obj.NumDataSubcarriers) ./ sqrt(obj.NumDataSubcarriers);
            if upper(obj.WaveformType) == "AFDM"
                chirpMatrix1 = diag(exp(-1j * 2 * pi * obj.ChirpParam1 * ((0:obj.NumDataSubcarriers-1).^2)));
                chirpMatrix2 = diag(exp(-1j * 2 * pi * obj.ChirpParam2 * ((0:obj.NumDataSubcarriers-1).^2)));
                % H_eff = A * H_time * A'
                % 其中 A = L2 * F * L1
                transformMatrix = chirpMatrix2 * dftMatrix * chirpMatrix1;
            else
                transformMatrix = dftMatrix;
            end

            
            effectiveChannelMatrix = transformMatrix * effectiveChannelTime * transformMatrix';

        end

        % --- 均衡器调度 ---
        function eqSignal = equalizer(obj, daftRx, effectiveChannelMatrix, noisePower)
            switch upper(obj.EqualizerType)
                case "MMSE"
                    eqSignal = runMmse(obj,daftRx, effectiveChannelMatrix, noisePower);
                case {"MRC","MRC-DFE"}
                    eqSignal = runWeightedMrcDfe(obj,daftRx, effectiveChannelMatrix, noisePower);
                otherwise
                    error('未知的均衡器类型: %s', obj.EqualizerType);
            end

        end

        % --- MMSE 均衡器 ---
        function estimatedSignal = runMmse(obj,receivedSignal, effectiveChannelMatrix, noisePower)
            
                % 只取有效载波对应的子矩阵
                effectiveChannelMatrixActive = effectiveChannelMatrix(:, obj.ActiveIndices); 
                % 计算 Gram 矩阵 (H' * H)
                gramMatrix = effectiveChannelMatrixActive' * effectiveChannelMatrixActive;   % Nr * Nr
                % 求解线性方程组 (H'H + N0*I) * x = H'y
                estimatedSignal = (gramMatrix + noisePower * eye(obj.NumActiveCarriers)) \ (effectiveChannelMatrixActive' * receivedSignal); 
        end

        % --- 加权 MRC-DFE 均衡器 ---
        function estimatedSignal = runWeightedMrcDfe(obj,receivedSignal, effectiveChannelMatrix, noisePower)
            
            % 功能: 基于加权最大比合并的判决反馈均衡器 (Weighted MRC-DFE) 信号检测
            % 输入:
            %   receivedSignal         : 去除了CPP的接收信号 (Nr x 1)
            %   effectiveChannelMatrix : 分数多普勒下的稀疏信道矩阵 (Nr x Nr)
            %   noisePower             : 噪声功率
            %   activeIndices          : 有效数据子载波的索引列表
            %   numIterations          : 迭代次数
            % 输出:
            %   estimatedSignal        : 估计的有用信号

            numDataSubcarriers = size(receivedSignal, 1);
            epsilon = 1e-5;
            
            % 预计算每列能量 (Column Energies) ||h_k||^2
            columnEnergies = full(sum(abs(effectiveChannelMatrix).^2, 1)).';
            
            % 稀疏矩阵索引加速 (缓存非零元素的行索引和值)
            columnRowIndices = cell(numDataSubcarriers, 1);
            columnValues = cell(numDataSubcarriers, 1);
            
            for k = 1 : numDataSubcarriers
                [rows, ~, values] = find(effectiveChannelMatrix(:, k));
                columnRowIndices{k} = rows;
                columnValues{k} = values;
            end
            
            currentEstimate = zeros(numDataSubcarriers, 1); % 初始化估计值 (ZP位置自动为0)
            previousEstimate = zeros(numDataSubcarriers, 1);
            residualSignal = receivedSignal; % 初始残差等于接收信号
            
            % --- 基于信道能量进行降序排序 ---
            activeEnergies = columnEnergies(obj.ActiveIndices);
            [~, sortIdx] = sort(activeEnergies, 'descend');
            orderedIndices = obj.ActiveIndices(sortIdx); % 能量大的子载波排在前面

            for n = 1 : obj.NumMaxIterations
                for k = orderedIndices
                    rows = columnRowIndices{k};
                    if isempty(rows)
                        continue; 
                    end
                    channelValues = columnValues{k};
                    
                    % 加权 MRC (Weighted MRC Step)
                    mrcOutput = sum(conj(channelValues) .* residualSignal(rows)) + columnEnergies(k) * previousEstimate(k);
                    
                    % LMMSE 更新
                    newEstimate = mrcOutput / (columnEnergies(k) + noisePower);
                    
                    % 计算变化量并更新残差
                    estimateChange = newEstimate - previousEstimate(k);
                    currentEstimate(k) = newEstimate; % 存储当前迭代的估计值
                    
                    % 如果变化足够大，则更新残差信号
                    if abs(estimateChange) > 1e-6 
                        residualSignal(rows) = residualSignal(rows) - channelValues * estimateChange;
                    end
                end

                % 检查本次迭代与上一次迭代的向量差异 (2-范数)
                iterationDiff = norm(currentEstimate - previousEstimate);
                
                % 如果变化量已经极小，说明算法已收敛，提前跳出循环
                if iterationDiff < epsilon
                    % 可以在这里取消注释下面这行来观察平均迭代次数
                    % fprintf('MRC-DFE 在第 %d 次迭代提前收敛\n', n);
                    break; 
                else
                    % 准备下一次迭代
                    previousEstimate = currentEstimate;
                end

            end
            
            % % 输出完整的向量（含ZP位置的0）
            % estimatedSignalFull = currentEstimate;

            % 直接输出裁剪好的有用数据
            estimatedSignal = currentEstimate(obj.ActiveIndices);
        end

        % --- 生成 DAFT 响应向量 ---
        function responseVector = buildDaftResponse(obj,currentDelayTap,currentDopplerShift)
            
            subcarrierIndices = 0 : obj.NumDataSubcarriers - 1;
            zeroBasedPilotIndex = obj.PilotIndex - 1;
            % q = obj.PilotIndex; % 导频位置索引 (在公式中通常记为 q)
            
            phaseTerm = (2 * pi / obj.NumDataSubcarriers) * (obj.NumDataSubcarriers * obj.ChirpParam1 * currentDelayTap^2 - ...
                            zeroBasedPilotIndex * currentDelayTap + obj.NumDataSubcarriers * obj.ChirpParam2 * (zeroBasedPilotIndex^2 - subcarrierIndices.^2));
            exponentialPhase = exp(1j * phaseTerm).';
            
            % 狄利克雷核
            thetaTerm = subcarrierIndices - zeroBasedPilotIndex + currentDopplerShift + 2 * obj.NumDataSubcarriers * obj.ChirpParam1 * currentDelayTap;
            
            numeratorTerm = exp(-1j * 2 * pi * thetaTerm) - 1;
            denominatorTerm = exp(-1j * 2 * pi / obj.NumDataSubcarriers * thetaTerm) - 1;
            
            dirichletKernel = zeros(obj.NumDataSubcarriers, 1);
            tol = 1e-6;
            
            % 处理分母为0的奇点
            singularityIndices = abs(denominatorTerm) < tol;
            dirichletKernel(~singularityIndices) = numeratorTerm(~singularityIndices) ./ denominatorTerm(~singularityIndices);
            dirichletKernel(singularityIndices) = obj.NumDataSubcarriers; 
            
            responseVector = (1 / obj.NumDataSubcarriers) * exponentialPhase .* dirichletKernel;
            
        end

        % --- 重建物理信道 ---
        function channelMatrixPhysical = rebuildLtvChannel(obj, pathParams)
            channelMatrixPhysical = sparse(obj.TotalSubcarriers, obj.TotalSubcarriers);
            numPaths = size(pathParams, 1);
            
            % 修正系数：从 N_data 域转换到 N_total 域
            dopplerCorrectionFactor = obj.TotalSubcarriers / obj.NumDataSubcarriers;

            for i = 1 : numPaths
                pathDelay = round(pathParams(i, 1));      % 时延 (l)
                pathDoppler = pathParams(i, 2);           % 归一化多普勒 (相对于 N_data)
                pathGain = pathParams(i, 3);              % 复增益 (h)
                
                % 转换多普勒基准
                pathDoppler = pathDoppler * dopplerCorrectionFactor; 
                
                % 构造线性移位矩阵 (Pi)
                rowIndices = (1 + pathDelay) : obj.TotalSubcarriers;
                colIndices = 1 : (obj.TotalSubcarriers - pathDelay);
                delayPermutationMatrix = sparse(rowIndices, colIndices, ones(1, length(rowIndices)), obj.TotalSubcarriers, obj.TotalSubcarriers);
                
                % 构造多普勒对角阵 (D)
                timeIndices = (0 : obj.TotalSubcarriers - 1).';
                dopplerDiagonal = exp(-1j * 2 * pi * pathDoppler * timeIndices / obj.TotalSubcarriers); 
                dopplerMatrix = spdiags(dopplerDiagonal, 0, obj.TotalSubcarriers, obj.TotalSubcarriers);
                
                % H = sum( h * D * P^l )
                channelMatrixPhysical = channelMatrixPhysical + pathGain * dopplerMatrix * delayPermutationMatrix;
            end

        end

    end
end
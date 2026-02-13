clc; clear all; close all;
%% 基础参数
    bitsPerSymbol = 2;                  % 每符号所含比特数
    modulationOrder = 2^bitsPerSymbol;  % 星座图点数
    numDatasSubcarriers = 2^8;          % 有效数据子载波数
    dataSnrRange = 0 : 5 : 20;          % 有效数据信噪比 (dB)
    pilotSnrRange = [25 30 35];         % 导频信噪比 (dB)
    numSNRPoints = length(dataSnrRange);

    numPaths = 3; % 多径数
    noisePower = 10.^(-dataSnrRange ./ 10);

    targetErrorCounts = 1e3;
    maxFrameCounts = 1e6;

    maxNormDoppler = 2; % 4GHz 载波下的 540 km/h
    % 给出固定的归一化时延，初始化归一化多普勒、信道增益
    pathDelays = [0 1 2];
    pathDopplers = zeros(numPaths);
    pathGains = zeros(numPaths);

    maxDelay = max(pathDelays);

    cppLength = maxDelay;

    carrierNum = numDatasSubcarriers + cppLength; % 所需子载波个数 Data + Overhead
    dataRange = (1 + cppLength) : carrierNum;

    % 设计c1、c2参数、设计帧结构 
    dopplerGuard = 2;    % 用于对抗分数多普勒，为不同径带来足够高的分离度
    chirpParam1 = (2 * (maxNormDoppler + dopplerGuard) + 1) / (2 * numDatasSubcarriers);           % 满足正交条件的c1，使多径各路径不发生交叠
    chirpParam2 = 1 / (numDatasSubcarriers^2 * 2 * pi);
    
    if (2 * (maxNormDoppler + dopplerGuard) * (maxDelay + 1)) + maxDelay > numDatasSubcarriers      %必须满足正交条件
        fprintf("子载波不满足正交！\n");
    end

    % ZP长度
    zeroPaddingLenth = (maxDelay + 1) * (2 * (maxNormDoppler + dopplerGuard) + 1) - 1;

    % 导频设计
    pilotIdx = 1;
    pilot = 10 + 0j;

    if pilotIdx + 2 * zeroPaddingLenth >= numDatasSubcarriers
        fprintf("dataNum过小,放不下双边间隔! \n");
    end
%% 蒙特卡洛仿真
    %% 多帧模拟
        bitErrorRates = zeros(1,numSNRPoints);    % 初始化BER数组
        %% SNR步进
        for i = 1 : numSNRPoints
            % 初始化
            totalErrs = 0;
            totalBits = 0;
            frameCounts = 0;

            fprintf("当前正在仿真 %d SNRdB ... \n", dataSnrRange(i));
            while(totalErrs <= targetErrorCounts && frameCounts < maxFrameCounts)
                % Jakes分布
                % P个随机相位theta
                theta = (rand(1, numPaths) * 2 * pi) - pi;
                pathDopplers = maxNormDoppler * cos(theta);
                
                % 瑞利信道
                pathGains = (randn(1, numPaths) + 1j * randn(1, numPaths)) / sqrt(2 * numPaths);

                maxDoppler = max(abs(pathDopplers)); 
                maxDoppler_alpha = round(maxDoppler);                 % 对多普勒取整
                maxDoppler_a = maxDoppler - maxDoppler_alpha;        % 多普勒的小数部分，范围[-0.5,0.5]
                    
                    % 有效数据个数
                    numActiveCarriers = numDatasSubcarriers - 2 * zeroPaddingLenth - 1;  % 减1是因为还要给导频留位置
                    activeIdx = (zeroPaddingLenth + 2) : (zeroPaddingLenth + 1 + numActiveCarriers); % 有效数据位置索引

                    if max(activeIdx) > (numDatasSubcarriers - zeroPaddingLenth)
                        error('索引分配错误');
                    end
                    % 构造LTV信道矩阵
                    channelMatrixPhysical = LtvChannel(carrierNum,pathDelays,pathDopplers,pathGains);

                    % 每生成一个符号，计数器+1
                    frameCounts = frameCounts + 1;
    
                %% 发送端
                    % 生成原始数据
                    originalData = randi([0 modulationOrder-1], numActiveCarriers, 1);
                    % 映射成星座点
                    qamSymbols = qammod(originalData, modulationOrder, 'UnitAveragePower', true);
                    % 组建完整的DAFT域符号向量 (加入ZP & pilot)
                    daftDomainFrame = zeros(numDatasSubcarriers,1);
                    daftDomainFrame(pilotIdx) = pilot;
                    daftDomainFrame(activeIdx) = qamSymbols;
                    % IDAFT(c1 = c2 = 0时IDAFT退化成为IDFT)
                    timeDomainFrame = IDAFT(daftDomainFrame, chirpParam1, chirpParam2);
                    % X = IDAFT(X_QAM, 0, 0);
                    % 生成CPP
                    cpp = timeDomainFrame(end - cppLength + 1 : end) .* exp(-1j * 2 * pi * chirpParam1 * (numDatasSubcarriers^2 + 2 * numDatasSubcarriers * (-cppLength:-1).'));
                    % CPP = X(end - CPP_lenth + 1 : end);
                    % 添加CPP / CP
                    txSignal = [cpp; timeDomainFrame];
                
                %% 叠加AWGN
                    noiseVec = sqrt(noisePower(i) / 2) * (randn(carrierNum,1) + 1j * randn(carrierNum,1));           % AWGN

                %% 接收端
                    % 接收到的时域符号
                    rxSignal = channelMatrixPhysical * txSignal + noiseVec;    % 接收信号 r = H * s + w
        
                    % 移除CPP
                    rxSignal = rxSignal(cppLength + 1 : end);

                    % DAFT (c1 = c2 = 0时DAFT退化成为DFT)
                    daftReceivedSignal = DAFT(rxSignal,chirpParam1,chirpParam2);

                    % 提取导频作信道估计
                    [channelMatrixEffective,est_params] = ChannelEstimator(daftReceivedSignal, pilot, pilotIdx, numDatasSubcarriers, cppLength, chirpParam1, chirpParam2, numPaths);


                    % 均衡器
                    % Y_hat_full = Weighted_MRC_DFE(Y,H_eff,N0(i),active_idx,10);
                    equalizedSignalFull = MmseEqualizer(daftReceivedSignal,channelMatrixEffective(:,activeIdx),noisePower(i));
                    equalizedSignalActive = equalizedSignalFull;

                    % 解映射
                    rxDemodulatedData = qamdemod(equalizedSignalActive, modulationOrder, 'UnitAveragePower', true);
        
                    % 统计帧错误
                    [frameErrors,~] = biterr(originalData,rxDemodulatedData);
                        
                    % 累加错误与统计总bit数
                    totalErrs = totalErrs + frameErrors;
                    totalBits = totalBits + numActiveCarriers * bitsPerSymbol;
                    % 实时显示进度，防止以为卡死
                    if mod(frameCounts, 1e3) == 0
                        fprintf('已跑 %d 帧, 收集错误数: %d\n',frameCounts, totalErrs);
                    end
            end
            bitErrorRates(i) = totalErrs / totalBits;
            fprintf('SNR %d dB 完成: 共跑 %d 帧, BER = %.2e\n', dataSnrRange(i), frameCounts, bitErrorRates(i));
        end

%% 绘图结果展示
    figure(1);
    semilogy(dataSnrRange,bitErrorRates,"-go","LineWidth",1.5);
    hold on;
    grid on;

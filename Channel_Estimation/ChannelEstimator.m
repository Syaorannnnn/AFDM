function [estimatedChannelMatrixEffective, finalEstimatedParams] = ChannelEstimator(receivedSignalDaft, pilotSymbol, pilotIndex, numDataSubcarriers, cppLength, chirpParam1, chirpParam2, numPaths)
    % channelEstimator : AFDM 分数多普勒信道估计与矩阵重建
    %
    % 输入:
    %   receivedSignalDaft : 接收到的 DAFT 域信号 (N_data x 1)
    %   pilotSymbol        : 发送的导频符号值 (标量)
    %   pilotIndex         : 导频在 DAFT 帧中的位置索引 (1 ~ N_data)
    %   numDataSubcarriers : 有效数据长度
    %   cppLength          : CPP 长度
    %   chirpParam1/2      : AFDM 参数
    %   numPaths           : 需要估计的路径数量
    %
    % 输出:
    %   estimatedChannelMatrixEffective : 重建后的等效信道矩阵
    %   finalEstimatedParams            : 估计参数 [Delay, Doppler, ComplexGain]

    % 搜索配置
    delaySearchRange = 0 : cppLength;              
    dopplerSearchRange = -2 : 2;
    dopplerStepSize = 1e-4;   % 精搜索步长            
    
    correlationGrid = zeros(length(delaySearchRange), length(dopplerSearchRange));
    tempEstimatedParams = zeros(numPaths, 3); 
    
    % 粗搜索
    for i = 1 : length(delaySearchRange)
        for j = 1 : length(dopplerSearchRange)
            currentDelay = delaySearchRange(i);
            currentIntegerDoppler = dopplerSearchRange(j);
            
            % pilotIndex - 1 用于将 MATLAB索引 转为 0-based 公式索引
            testResponseVector = BuildDaftResponse(numDataSubcarriers, chirpParam1, chirpParam2, currentDelay, currentIntegerDoppler, pilotIndex - 1);
            
            correlationGrid(i, j) = abs(testResponseVector' * receivedSignalDaft);
        end
    end
    
    tempCorrelationGrid = correlationGrid;
    
    % 迭代提取路径
    for p = 1 : numPaths
        [~, maxLinearIndex] = max(tempCorrelationGrid(:));
        [rowIndex, colIndex] = ind2sub(size(tempCorrelationGrid), maxLinearIndex);
        
        estimatedDelay = delaySearchRange(rowIndex);
        estimatedIntegerDoppler = dopplerSearchRange(colIndex);
        
        % 精搜索
        bestCorrelation = -inf;
        bestFractionalDoppler = 0;
        
        for fractionalDoppler = -0.5 : dopplerStepSize : 0.5
            currentDopplerTest = estimatedIntegerDoppler + fractionalDoppler;
            fineResponseVector = BuildDaftResponse(numDataSubcarriers, chirpParam1, chirpParam2, estimatedDelay, currentDopplerTest, pilotIndex - 1);
            
            currentCorrelation = abs(fineResponseVector' * receivedSignalDaft);
            if currentCorrelation > bestCorrelation
                bestCorrelation = currentCorrelation;
                bestFractionalDoppler = fractionalDoppler;
            end
        end
        estimatedDoppler = estimatedIntegerDoppler + bestFractionalDoppler;
        
        % 暂存参数 (Delay, Doppler, 占位Gain)
        tempEstimatedParams(p, :) = [estimatedDelay, estimatedDoppler, 0];
        
        % 屏蔽邻域
        colStart = max(1, colIndex - 2); 
        colEnd   = min(size(tempCorrelationGrid, 2), colIndex + 2);
        tempCorrelationGrid(rowIndex, colStart : colEnd) = -inf; 
    end
    
    % LS 求解增益
    basisMatrix = zeros(numDataSubcarriers, numPaths);
    for p = 1 : numPaths
        delayP = tempEstimatedParams(p, 1);
        dopplerP = tempEstimatedParams(p, 2);
        basisMatrix(:, p) = BuildDaftResponse(numDataSubcarriers, chirpParam1, chirpParam2, delayP, dopplerP, pilotIndex - 1);
    end
    
    systemMatrix = basisMatrix * pilotSymbol;
    optimizedGains = pinv(systemMatrix) * receivedSignalDaft;
    
    tempEstimatedParams(:, 3) = optimizedGains;
    
    % 排序与相位去旋转
    % 按时延排序，方便对比
    [~, sortIndices] = sort(tempEstimatedParams(:, 1)); 
    sortedEstimatedParams = tempEstimatedParams(sortIndices, :);
    
    % 去旋转,补偿 CPP 移除造成的相位偏移
    derotationPhasor = exp(1j * 2 * pi * sortedEstimatedParams(:, 2) * (cppLength / numDataSubcarriers));
    
    finalEstimatedParams = sortedEstimatedParams;
    finalEstimatedParams(:, 3) = sortedEstimatedParams(:, 3) .* derotationPhasor;
    
    % 重建信道
    totalSubcarriers = numDataSubcarriers + cppLength;
    
    % 重建物理信道 (LTV)
    tempPhysicalChannelMatrix = RebuildLtvChannel(totalSubcarriers, numDataSubcarriers, finalEstimatedParams);
    
    % 生成等效信道
    estimatedChannelMatrixEffective = GenerateEffectiveChannelMatrix(tempPhysicalChannelMatrix, cppLength, chirpParam1, chirpParam2);

end
function [estimatedSignalFull] = WeightedMrcDfe(receivedSignal, effectiveChannelMatrix, noisePower, activeIndices, numIterations)
    
    % 功能: 基于加权最大比合并的判决反馈均衡器 (Weighted MRC-DFE) 信号检测
    % 输入:
    %   receivedSignal         : 去除了CPP的接收信号 (Nr x 1)
    %   effectiveChannelMatrix : 分数多普勒下的稀疏信道矩阵 (Nr x Nr)
    %   noisePower             : 噪声功率
    %   activeIndices          : 有效数据子载波的索引列表
    %   numIterations          : 迭代次数
    % 输出:
    %   estimatedSignalFull    : 完整的估计信号向量 (含ZP部分的0)

    numDataSubcarriers = size(receivedSignal, 1);
    
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
    
    for n = 1 : numIterations
        % === 核心修改：只遍历有效数据索引 ===
        for k = activeIndices 
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
        % 准备下一次迭代
        previousEstimate = currentEstimate;
    end
    
    % 输出完整的向量（含ZP位置的0）
    estimatedSignalFull = currentEstimate;
end
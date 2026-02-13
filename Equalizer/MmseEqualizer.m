function [estimatedSignal] = MmseEqualizer(receivedSignal, effectiveChannelMatrix, noisePower)
    
    % 功能： MMSE 线性均衡器 复杂度O(N^3)
    % 输入:
    %   receivedSignal         : 接收的已经去除CPP的信号 (Nr x 1)
    %   effectiveChannelMatrix : 等效信道矩阵 (Nr x Nr)
    %   noisePower             : 噪声功率
    % 输出:
    %   estimatedSignal        : 均衡后的信号估计值

    [~, numSubcarriers] = size(effectiveChannelMatrix);

    % 计算 Gram 矩阵 (H' * H)
    gramMatrix = effectiveChannelMatrix' * effectiveChannelMatrix;   % Nr * Nr
 
    % 求解线性方程组 (H'H + N0*I) * x = H'y
    estimatedSignal = (gramMatrix + noisePower * eye(numSubcarriers)) \ (effectiveChannelMatrix' * receivedSignal); 

    % fprintf("%.2e\n", cond(gramMatrix + noisePower * eye(numSubcarriers)));
end
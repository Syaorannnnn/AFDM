function [channelMatrixEffective] = GenerateEffectiveChannelMatrix(channelMatrixPhysical, cppLength, chirpParam1, chirpParam2)

% 功能:   等效信道矩阵生成
% 输入：
%   channelMatrixPhysical ：双色散物理信道矩阵 (H)
%   cppLength             ：CPP长度
%   chirpParam1           ：AFDM 参数 c1
%   chirpParam2           ：AFDM 参数 c2
% 输出：
%   channelMatrixEffective: AFDM等效信道矩阵 (H_eff)

    totalSubcarriers = size(channelMatrixPhysical, 1);
    numDataSubcarriers = totalSubcarriers - cppLength;
    dataIndices = (1 + cppLength) : totalSubcarriers;
    
    % 构造 CPP 添加矩阵 (M)
    cppInsertionMatrix = zeros(totalSubcarriers, numDataSubcarriers);
    
    % 数据部分直接映射
    cppInsertionMatrix(cppLength + 1 : end, :) = eye(numDataSubcarriers);
    
    % CPP 部分 (基于数据尾部生成)
    gammaVector = exp(-1j * 2 * pi * chirpParam1 * (numDataSubcarriers^2 + 2 * numDataSubcarriers * (-cppLength:-1).'));
    cppInsertionMatrix(1 : cppLength, (numDataSubcarriers - cppLength + 1) : numDataSubcarriers) = diag(gammaVector);

    % 计算等效时域矩阵 (N_data x N_data)
    % 截取接收信号的数据部分 = H(dataIndices, :) * M * S0
    effectiveChannelTime = channelMatrixPhysical(dataIndices, :) * cppInsertionMatrix;

    % 变换域矩阵
    chirpMatrix1 = diag(exp(-1j * 2 * pi * chirpParam1 * ((0:numDataSubcarriers-1).^2)));
    chirpMatrix2 = diag(exp(-1j * 2 * pi * chirpParam2 * ((0:numDataSubcarriers-1).^2)));
    dftMatrix  = dftmtx(numDataSubcarriers) ./ sqrt(numDataSubcarriers);

    % H_eff = A * H_time * A'
    % 其中 A = L2 * F * L1
    transformMatrix = chirpMatrix2 * dftMatrix * chirpMatrix1;
    
    channelMatrixEffective = transformMatrix * effectiveChannelTime * transformMatrix';

end

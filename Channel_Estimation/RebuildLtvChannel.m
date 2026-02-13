function channelMatrixPhysical = RebuildLtvChannel(totalSubcarriers, numDataSubcarriers, pathParams)
    % rebuildLtvChannel : 重建物理信道 (线性移位模型)
    % pathParams: [delay, doppler, gain]
    
    channelMatrixPhysical = sparse(totalSubcarriers, totalSubcarriers);
    numPaths = size(pathParams, 1);
    
    % 修正系数：从 N_data 域转换到 N_total 域
    dopplerCorrectionFactor = totalSubcarriers / numDataSubcarriers;

    for i = 1 : numPaths
        pathDelay = round(pathParams(i, 1));      % 时延 (l)
        pathDoppler = pathParams(i, 2);           % 归一化多普勒 (相对于 N_data)
        pathGain = pathParams(i, 3);              % 复增益 (h)
        
        % 转换多普勒基准
        pathDoppler = pathDoppler * dopplerCorrectionFactor; 
        
        % 构造线性移位矩阵 (Pi)
        rowIndices = (1 + pathDelay) : totalSubcarriers;
        colIndices = 1 : (totalSubcarriers - pathDelay);
        delayPermutationMatrix = sparse(rowIndices, colIndices, ones(1, length(rowIndices)), totalSubcarriers, totalSubcarriers);
        
        % 构造多普勒对角阵 (D)
        timeIndices = (0 : totalSubcarriers - 1).';
        dopplerDiagonal = exp(-1j * 2 * pi * pathDoppler * timeIndices / totalSubcarriers); 
        dopplerMatrix = spdiags(dopplerDiagonal, 0, totalSubcarriers, totalSubcarriers);
        
        % H = sum( h * D * P^l )
        channelMatrixPhysical = channelMatrixPhysical + pathGain * dopplerMatrix * delayPermutationMatrix;
    end

end
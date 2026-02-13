function channelMatrixPhysical = LtvChannel(totalSubcarriers, pathDelays, pathDopplers, pathGains)
    
    % 功能：生成线性时变信道矩阵 (LTV)
    % 输入：
    %   totalSubcarriers : 总子载波数 (含CPP)
    %   pathDelays       : 路径时延数组
    %   pathDopplers     : 路径多普勒数组
    %   pathGains        : 路径增益数组

    numPaths = length(pathDelays);

    % 构造循环移位矩阵基 (用于模拟时延)
    % 虽然叫 LTV，但这里的 Pi 构造是基于循环移位的 (Toeplitz)
    baseVector = [zeros(1, totalSubcarriers - 1) 1];
    delayPermutationMatrix = toeplitz([baseVector(1) fliplr(baseVector(2:end))], baseVector);

    channelMatrixPhysical = zeros(totalSubcarriers, totalSubcarriers);

    for i = 1 : numPaths
        pathGain = pathGains(i);
        pathDelay = pathDelays(i);     % 时延
        pathDoppler = pathDopplers(i); % 多普勒
        
        % 多普勒相位矩阵 (对角阵)
        dopplerMatrix = diag(exp(-1j * 2 * pi * pathDoppler * (0 : totalSubcarriers - 1) / totalSubcarriers));

        % H = sum( h_i * D_i * P^l_i )
        channelMatrixPhysical = channelMatrixPhysical + pathGain * dopplerMatrix * (delayPermutationMatrix^pathDelay);
    end

end
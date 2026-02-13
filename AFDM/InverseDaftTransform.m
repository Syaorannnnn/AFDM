function timeDomainSignal = InverseDaftTransform(daftDomainSignal, chirpParam1, chirpParam2)
        
    % 功能： 离散仿射傅里叶反变换
    % 输入：
    %   daftDomainSignal ：DAFT域原始数据符号 (X)
    %   chirpParam1      ：线性调制信号的第一个参数 (c1)
    %   chirpParam2      ：线性调制信号的第二个参数 (c2)
    % 输出：
    %   timeDomainSignal ：时域已调信号 (S)
    % 

    signalLength = size(daftDomainSignal, 1);

    chirpMatrix1 = diag(exp(-1j * 2 * pi * chirpParam1 * ((0:signalLength-1).^2)));
    chirpMatrix2 = diag(exp(-1j * 2 * pi * chirpParam2 * ((0:signalLength-1).^2)));
    dftMatrix = dftmtx(signalLength) ./ sqrt(signalLength);

    transformMatrix = chirpMatrix2 * dftMatrix * chirpMatrix1;
    % IDAFT 是 DAFT 的逆（对于酉矩阵是共轭转置）
    timeDomainSignal = transformMatrix' * daftDomainSignal;

end

function demodulatedSignal = DaftTransform(receivedSignal, chirpParam1, chirpParam2)
    
    % 功能： 离散仿射傅里叶变换
    % 输入：
    %   receivedSignal ：接收端接收到的符号
    %   chirpParam1    ：线性调制信号的第一个参数 (c1)
    %   chirpParam2    ：线性调制信号的第二个参数 (c2)
    % 输出：
    %   demodulatedSignal：解调结果 (Y)
    %   Y = A*r = L2*F*L1*r

    signalLength = size(receivedSignal, 1);

    chirpMatrix1 = diag(exp(-1j * 2 * pi * chirpParam1 * ((0:signalLength-1).^2)));
    chirpMatrix2 = diag(exp(-1j * 2 * pi * chirpParam2 * ((0:signalLength-1).^2)));
    dftMatrix = dftmtx(signalLength) ./ sqrt(signalLength);

    transformMatrix = chirpMatrix2 * dftMatrix * chirpMatrix1;
    demodulatedSignal = transformMatrix * receivedSignal;

end

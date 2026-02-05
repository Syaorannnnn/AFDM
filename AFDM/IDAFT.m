function S = IDAFT(X, c1, c2)
        
    % 功能： 离散仿射傅里叶反变换
    % 输入：
    %   X   ：原始数据符号
    %   c1  ：线性调制信号的第一个参数
    %   c2  ：线性调制信号的第一个参数
    % 输出：
    %   S   ：已调信号 
    % 

    N = size(X,1);

    L1 = diag(exp(-1j * 2 * pi * c1 * ((0:N-1).^2)));
    L2 = diag(exp(-1j * 2 * pi * c2 * ((0:N-1).^2)));
    F =  dftmtx(N) ./ sqrt(N);

    A = L2 * F * L1;
    S = A' * X;

end

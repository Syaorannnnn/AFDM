function Y = DAFT(R, c1, c2)
    
% 功能： 离散仿射傅里叶变换
% 输入：
%   R   ：接收端接收到的符号
%   c1  ：线性调制信号的第一个参数
%   c2  ：线性调制信号的第一个参数
% 输出：
%   Y：解调结果 
%   Y = A*r = L2*F*L1*r

N = size(R,1);

L1 = diag(exp(-1j * 2 * pi * c1 * ((0:N-1).^2)));
L2 = diag(exp(-1j * 2 * pi * c2 * ((0:N-1).^2)));
F =  dftmtx(N) ./ sqrt(N);

A = L2 * F * L1;
Y = A * R;

end

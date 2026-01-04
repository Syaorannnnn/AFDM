function Heff = AFDM_Channel(H, c1, c2)

% 功能:   AFDM等效信道矩阵生成
% 输入：
%   H   ：双色散信道
%   c1  ：线性调制信号的第一个参数
%   c2  ：线性调制信号的第一个参数
% 输出：
%   Heff: AFDM等效信道矩阵

N = size(H,1);
n = (0:N-1).';

L1 = diag(exp(-1j * 2 * pi * c1 * (n.^2)));
L2 = diag(exp(-1j * 2 * pi * c2 * (n.^2)));
F  = dftmtx(N) ./ sqrt(N);

% Heff = L2 * F * L1 * H * L1' * F' * L2'
Heff = L2 * F * L1 * H * L1' * F' * L2';

end

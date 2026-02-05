function [H_eff] = Gen_Eff_Channel(H,CPP_lenth, c1, c2)

% 功能:   等效信道矩阵生成
% 输入：
%   H   ：双色散信道
%   CPP_lenth: CPP长度
%   c1  ：线性调制信号的第一个参数
%   c2  ：线性调制信号的第一个参数
% 输出：
%   H_eff: AFDM等效信道矩阵

    N = size(H,1);
    N_data = N - CPP_lenth;
    data_range = (1 + CPP_lenth) : N;
    M = zeros(N, N_data);
    % 数据部分直接映射
    M(CPP_lenth + 1 : end, :) = eye(N_data);
    % CPP / CP 部分
    Gamma = exp(-1j * 2 * pi * c1 * (N_data^2 + 2 * N_data * (-CPP_lenth:-1).'));
    M(1 : CPP_lenth, (N_data - CPP_lenth + 1) : N_data) = diag(Gamma);

    % 计算等效时域矩阵 (N_data x N_data)
    % 物理接收 = H * M * S0
    % 截取数据部分 = H(data_range, :) * M * S0
    % H_time_eff = H(data_range, :) * M
    H_eff = H(data_range, :) * M;

    L1 = diag(exp(-1j * 2 * pi * c1 * ((0:N_data-1).^2)));
    L2 = diag(exp(-1j * 2 * pi * c2 * ((0:N_data-1).^2)));
    F  = dftmtx(N_data) ./ sqrt(N_data);

    H_eff = L2 * F * L1 * H_eff * L1' * F' * L2';

end

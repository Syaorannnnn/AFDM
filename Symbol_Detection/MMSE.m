function [X_hat] = MMSE(R, Heff, N0)
    
    % 功能： MMSE 复杂度O(N^3)
    % 输入:
    %   R           : 接收的已经去除CPP的信号 (Nr x 1)
    %   Heff        : 等效信道矩阵 (Nr x Nr)
    %   N0          : 噪声功率

    [~,Nr] = size(Heff);

    % 厄密特矩阵
    Q = Heff' * Heff;   % Nr * Nr
 
    % 求解方程组
    X_hat = (Q + N0 * eye(Nr)) \ (Heff' * R); 

    % fprintf("%.2e\n",cond(Q + N0*eye(Nr)));
end
function [X_hat] = Iterative_MMSE(y, Heff_t, N0, n_iter)
    % 功能： 使用 Gauss-Seidel 迭代法实现 MMSE 检测
    % 输入:
    %   y           : 接收信号 (Nx1)
    %   Heff_t      : 分数多普勒下的稀疏信道矩阵 (NxN)
    %   N0          : 噪声功率
    %   n_iter      : 迭代次数

    [~, K] = size(Heff_t);

    % 厄密特矩阵
    Q = Heff_t' * Heff_t;
    b = Heff_t' * y;

    D = diag(Q);        % Q的主对角线
    L1 = tril(Q, -1);   % Q不含主对角线的下三角部分
    L2 = triu(Q, 1);    % Q不含主对角线的上三角部分
    
    D = D + N0;       % 总对角线 = 信号能量 + 噪声
    A = L1 + diag(D); 
    
    X_hat = zeros(K, 1);

    for n = 1 : n_iter
        % 计算残差
        RHS = b - L2 * X_hat;
        
        % 求解方程组
        X_hat = A \ RHS; 
    end

end
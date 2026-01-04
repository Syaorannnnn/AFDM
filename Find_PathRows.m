function pathRows = Find_PathRows(H, K)
    % 简单的全能量搜索，每列找最强的 K 个点
    [~, N_sym] = size(H);
    pathRows = zeros(K, N_sym);
    for i = 1 : N_sym
        [~, idx] = sort(abs(H(:, i)), 'descend');
        pathRows(:, i) = idx(1:K);
    end
end

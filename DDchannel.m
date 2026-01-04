function [H] = DDchannel(N, pathDelays, pathDopplers, pathGains,samplingFreq,c1)

% 双色散信道模型
% 输入：
%   N               : 符号长度
%   pathDelays      : 时延
%   pathDopplers    : 多普勒频偏
%   pathGains       : 复增益
% 输出：
%   H               : 信道矩阵

%% 归一化时延/多普勒
P = length(pathDelays);
% normDelays = round(pathDelays * samplingFreq);      % 归一化延迟为整数 li
% normDopplers = pathDopplers / (N * samplingFreq);     % 归一化的多普勒频移(数字频率) vi = N * fi = αi + ai
normDelays = pathDelays;            % 归一化延迟为整数 li
normDopplers = pathDopplers;        % 归一化的多普勒频移(数字频率) vi = N * fi = αi + ai
%fprintf("%s",normDopplers);
%% 矩阵构造
% PI
Pi = [zeros(1,N - 1) 1];
Pi = toeplitz([Pi(1) fliplr(Pi(2:end))], Pi);
% 
% %Delta_fi = W^normDopplers(i)
% W = diag(exp(-1j * 2 * pi * (0 : N-1) / N));
% 
% % Chirp 导频矩阵 Gamma_CPP
% H_p = zeros(N,N,P);
% CPP_p_temp = zeros(N,N,P);
% for p = 1 : P
%     CPP_p = diag( [exp(-1j * 2 * pi * c1 * ((N^2) + 2 * N * ((0:(normDelays(p) - 1)) - normDelays(p)))) , ones(1, N - normDelays(p))]);
%     % 第p路信道矩阵
%     H_p(:,:,p) = pathGains(p) * CPP_p * W^normDopplers(p) * Pi^normDelays(p);
%     CPP_p_temp(:,:,p) = CPP_p;
% end

H = zeros(N,N);
for i = 1 : P
    h_i = pathGains(i);
    l_i = normDelays(i);
    f_i = normDopplers(i);
    D_i = diag(exp(-1j * 2 * pi * f_i * (0 : N - 1) / N));
    for n = 0 : N-1
        if n < l_i
            temp(n+1) = exp(-1j * 2 * pi * c1 * (N^2 - 2 * N * (l_i - n)));
        else
            temp(n+1) = 1;
        end
    end
    G_i = diag(temp);   % equation (26) in [R1]
    H = H + h_i * G_i * D_i * Pi^l_i;
end

% % 完整信道矩阵 P路 (时延 * 多普勒 * 增益) 求和
% figure(4)
% subplot(1,3,1)
% imagesc(abs(H_p(:,:,1)))
% subplot(1,3,2)
% imagesc(abs(H_p(:,:,2)))
% subplot(1,3,3)
% imagesc(abs(H_p(:,:,3)))
% 
% H = sum(H_p,3);
% CPP = sum(CPP_p_temp,3);

end
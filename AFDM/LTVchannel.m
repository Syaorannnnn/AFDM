function H = LTVchannel(N,pathDelays,pathDopplers,pathGains)

P = length(pathDelays);
normDelays = pathDelays;            % 假设输入已经是整数采样点
normDopplers = pathDopplers;        % 归一化多普勒

Pi = [zeros(1,N - 1) 1];
Pi = toeplitz([Pi(1) fliplr(Pi(2:end))], Pi);

H = zeros(N,N);

for i = 1 : P
    h_i = pathGains(i);
    l_i = normDelays(i);     % 时延
    f_i = normDopplers(i);   % 多普勒
    
    D_i = diag(exp(-1j * 2 * pi * f_i * (0 : N - 1) / N));

    H = H + h_i * D_i * Pi^l_i;
end

end
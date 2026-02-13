clc; clear; close all;

%% 实例化对象
afdm1 = AfdmSystem(); 
afdm2 = AfdmSystem();
%% 仿真配置
snrRange = 0 : 5 : 20;
berResults1 = zeros(size(snrRange));
berResults2 = zeros(size(snrRange));
afdm1.EqualizerType = "MMSE";

afdm2.EqualizerType = "MRC";
afdm2.NumIterations = 3;
%% 运行仿真
for i = 1 : length(snrRange)
    berResults1(i) = afdm1.runMonteCarlo(snrRange(i), 1e3, 1e6);
end

for i = 1 : length(snrRange)
    berResults2(i) = afdm2.runMonteCarlo(snrRange(i), 1e3, 1e6);
end
%%  绘图
figure;
semilogy(snrRange, berResults1, '-go',snrRange, berResults2,'-b^','LineWidth', 1.5);
grid on;
xlabel('SNR_d (dB)');
ylabel('BER');
title('AFDM BER Performance');
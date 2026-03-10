classdef AfdmTransforms
    % AfdmTransforms: AFDM 核心变换工具类
    % 近期改动：
    %   1. generateEffectiveChannelMatrix 中消除 dftmtx 和 diag:
    %      - diag(vec) * M  → vec .* M      (O(N²) 代替 O(N³))
    %      - M * diag(vec)  → M .* vec.'     (O(N²) 代替 O(N³))
    %      - dftmtx(N) * M  → fft(M)/√N      (O(NlogN) 代替 O(N²))
    %   2. daft/idaft 已使用向量化 chirp + fft (无需改动)

    methods (Static)
        %% DAFT / IDAFT

        function demodulatedSignal = daft(inputSignal, chirpParam1, chirpParam2)
            signalLength = size(inputSignal, 1);
            timeVector = (0:signalLength - 1).';
            chirpVector1 = exp(-1j * 2 * pi * chirpParam1 * (timeVector .^ 2));
            chirpVector2 = exp(-1j * 2 * pi * chirpParam2 * (timeVector .^ 2));

            % 支持 inputSignal 为多列矩阵 (隐式扩展)
            tempSignal = inputSignal .* chirpVector1;
            tempSignal = fft(tempSignal) ./ sqrt(signalLength);
            demodulatedSignal = tempSignal .* chirpVector2;
        end

        function timeDomainSignal = idaft(inputSignal, chirpParam1, chirpParam2)
            signalLength = size(inputSignal, 1);
            timeVector = (0:signalLength - 1).';
            chirpVector1 = exp(-1j * 2 * pi * chirpParam1 * (timeVector .^ 2));
            chirpVector2 = exp(-1j * 2 * pi * chirpParam2 * (timeVector .^ 2));

            tempSignal = inputSignal .* conj(chirpVector2);
            tempSignal = ifft(tempSignal) .* sqrt(signalLength);
            timeDomainSignal = tempSignal .* conj(chirpVector1);
        end

        %% DFT / IDFT

        function outputSignal = dft(inputSignal)
            signalLength = size(inputSignal, 1);
            outputSignal = fft(inputSignal) ./ sqrt(signalLength);
        end

        function outputSignal = idft(inputSignal)
            signalLength = size(inputSignal, 1);
            outputSignal = ifft(inputSignal) .* sqrt(signalLength);
        end

        %% 等效信道矩阵 (性能优化版)

        function effectiveChannelMatrix = generateEffectiveChannelMatrix(physicalChannelMatrix, configParams)
            % generateEffectiveChannelMatrix - 优化版
            %
            % 计算: H_eff = T · H_time · T'
            %   其中 T = chirpMatrix2 · dftMatrix · chirpMatrix1
            %
            % 优化策略:
            %   原版使用 dftmtx(N) 构造 N×N 满矩阵 + diag() 构造对角满矩阵,
            %   三次 N×N 矩阵乘法, 总计 O(N³) 且内存开销 O(N²).
            %
            %   优化版将每次乘法拆解为逐元素操作 + FFT/IFFT:
            %   左乘 T·M:
            %     step1: chirpVec1 .* M         (逐行乘, 替代 chirpMatrix1 * M)
            %     step2: fft(step1) / √N        (替代 dftMatrix * step1)
            %     step3: chirpVec2 .* step2      (替代 chirpMatrix2 * step2)
            %
            %   右乘 M·T':
            %     T' = chirpMatrix1' · dftMatrix' · chirpMatrix2'
            %     step1: M .* conj(chirpVec1).'  (右乘 chirpMatrix1' = 每列乘标量)
            %     step2: ifft(step1.').' · √N    (右乘 dftMatrix' = ifft 作用于行)
            %     step3: step2 .* conj(chirpVec2).' (右乘 chirpMatrix2')

            N = configParams.NumDataSubcarriers;
            prefixLength = configParams.PrefixLength;
            totalSubcarriers = configParams.TotalSubcarriers;

            dataIndices = (1 + prefixLength):totalSubcarriers;

            % --- 插入矩阵 (含 CPP 结构, 此处无法避免显式构造) ---
            insertionMatrix = zeros(totalSubcarriers, N);
            insertionMatrix(prefixLength + 1:end, :) = eye(N);

            if upper(configParams.WaveformType) == "AFDM"
                gammaVector = exp(-1j * 2 * pi * configParams.ChirpParam1 * ...
                    (N ^ 2 + 2 * N * (-prefixLength:-1).'));
            else
                gammaVector = ones(prefixLength, 1);
            end

            insertionMatrix(1:prefixLength, (N - prefixLength + 1):N) = diag(gammaVector);

            % --- 时域等效信道 (与原版一致) ---
            effectiveTimeChannel = physicalChannelMatrix(dataIndices, :) * insertionMatrix;
            % effectiveTimeChannel 是 N × N 矩阵

            % --- 核心优化: 替代 transformMatrix 的构造与乘法 ---
            if upper(configParams.WaveformType) == "AFDM"
                timeVector = (0:N - 1).';
                chirpVec1 = exp(-1j * 2 * pi * configParams.ChirpParam1 * (timeVector .^ 2)); % N×1
                chirpVec2 = exp(-1j * 2 * pi * configParams.ChirpParam2 * (timeVector .^ 2)); % N×1

                % 左乘 transformMatrix * effectiveTimeChannel
                temp = chirpVec1 .* effectiveTimeChannel; % diag(chirp1) * M
                temp = fft(temp) ./ sqrt(N); % dftMatrix * M
                temp = chirpVec2 .* temp; % diag(chirp2) * M

                % 右乘 temp * transformMatrix'
                % T' = chirpMatrix1^H · dftMatrix^H · chirpMatrix2^H
                % 右乘时按从左到右的顺序: 先 chirpMatrix1^H, 再 dftMatrix^H, 最后 chirpMatrix2^H
                temp = temp .* conj(chirpVec1).'; % 右乘 chirpMatrix1^H
                temp = ifft(temp.').'; % 右乘 dftMatrix^H (ifft 逐列→转置)
                temp = temp .* sqrt(N); % ifft 归一化补偿
                effectiveChannelMatrix = temp .* conj(chirpVec2).'; % 右乘 chirpMatrix2^H
            else
                % OFDM 模式: T = dftMatrix, 直接用 fft/ifft
                temp = fft(effectiveTimeChannel) ./ sqrt(N); % 左乘 dftMatrix
                effectiveChannelMatrix = (ifft(temp.') .* sqrt(N)).'; % 右乘 dftMatrix'
            end

        end

    end

end

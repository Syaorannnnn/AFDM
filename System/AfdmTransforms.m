classdef AfdmTransforms
    % AfdmTransforms: 包含 AFDM 系统所有底层数学变换的静态方法库

    methods (Static)

        % --- DAFT 正变换 ---
        function demodulatedSignal = daft(signal, c1, c2)
            signalLength = size(signal, 1);
            chirpMatrix1 = diag(exp(-1j * 2 * pi * c1 * ((0:signalLength - 1) .^ 2)));
            chirpMatrix2 = diag(exp(-1j * 2 * pi * c2 * ((0:signalLength - 1) .^ 2)));
            dftMatrix = dftmtx(signalLength) ./ sqrt(signalLength);

            transformMatrix = chirpMatrix2 * dftMatrix * chirpMatrix1;
            demodulatedSignal = transformMatrix * signal;
        end

        % --- IDAFT 逆变换 ---
        function timeDomainSignal = idaft(signal, c1, c2)
            signalLength = size(signal, 1);
            chirpMatrix1 = diag(exp(-1j * 2 * pi * c1 * ((0:signalLength - 1) .^ 2)));
            chirpMatrix2 = diag(exp(-1j * 2 * pi * c2 * ((0:signalLength - 1) .^ 2)));
            dftMatrix = dftmtx(signalLength) ./ sqrt(signalLength);

            transformMatrix = chirpMatrix2 * dftMatrix * chirpMatrix1;
            timeDomainSignal = transformMatrix' * signal;
        end

        % ---  DFT / IDFT ---
        function outSignal = dft(signal)
            signalLength = size(signal, 1);
            outSignal = (dftmtx(signalLength) ./ sqrt(signalLength)) * signal;
        end

        function outSignal = idft(signal)
            signalLength = size(signal, 1);
            outSignal = (dftmtx(signalLength) ./ sqrt(signalLength))' * signal;
        end

        % --- 等效信道矩阵生成 ---
        function effectiveChannelMatrix = generateEffectiveChannelMatrix(physicalChannelMatrix, config)
            numDataSubcarriers = config.NumDataSubcarriers;
            prefixLength = config.PrefixLength;
            totalSubcarriers = config.TotalSubcarriers;

            dataIndices = (1 + prefixLength):totalSubcarriers;

            % 构造 CPP/CP 添加矩阵 (M)
            insertionMatrix = zeros(totalSubcarriers, numDataSubcarriers);
            insertionMatrix(prefixLength + 1:end, :) = eye(numDataSubcarriers);

            if upper(config.WaveformType) == "AFDM"
                % CPP 的 gamma 补偿相位
                gammaVector = exp(-1j * 2 * pi * config.ChirpParam1 * (numDataSubcarriers ^ 2 + 2 * numDataSubcarriers * (-prefixLength:-1).'));
            else
                % OFDM 的普通 CP
                gammaVector = ones(prefixLength, 1);
            end

            insertionMatrix(1:prefixLength, (numDataSubcarriers - prefixLength + 1):numDataSubcarriers) = diag(gammaVector);

            % 等效时域矩阵
            effectiveTimeChannel = physicalChannelMatrix(dataIndices, :) * insertionMatrix;

            % 变换到 DAFT/DFT 域
            dftMatrix = dftmtx(numDataSubcarriers) ./ sqrt(numDataSubcarriers);

            if upper(config.WaveformType) == "AFDM"
                chirp1 = diag(exp(-1j * 2 * pi * config.ChirpParam1 * ((0:numDataSubcarriers - 1) .^ 2)));
                chirp2 = diag(exp(-1j * 2 * pi * config.ChirpParam2 * ((0:numDataSubcarriers - 1) .^ 2)));
                transformMatrix = chirp2 * dftMatrix * chirp1;
            else
                transformMatrix = dftMatrix;
            end

            rawEffectiveChannelMatrix = transformMatrix * effectiveTimeChannel * transformMatrix';

            % 检查 config 中是否配置了成形窗，若有则应用双边对角阵相乘
            if isprop(config, 'pulseShapingWindow') && ~isempty(config.PulseShapingWindow)
                windowDiagonal = diag(config.PulseShapingWindow);
                effectiveChannelMatrix = windowDiagonal * rawEffectiveChannelMatrix * windowDiagonal;
            else
                effectiveChannelMatrix = rawEffectiveChannelMatrix;
            end

        end

    end

end

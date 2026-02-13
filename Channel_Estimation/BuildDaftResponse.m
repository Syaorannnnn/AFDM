function responseVector = BuildDaftResponse(numSubcarriers, chirpParam1, chirpParam2, delayTap, dopplerShift, pilotIndex)
    % buildDaftResponse: 生成 AFDM 理论响应向量
    
    subcarrierIndices = 0 : numSubcarriers-1; 
    q = pilotIndex; % 导频位置索引 (在公式中通常记为 q)
    
    phaseTerm = (2*pi/numSubcarriers) * (numSubcarriers*chirpParam1*delayTap^2 - q*delayTap + numSubcarriers*chirpParam2*(q^2 - subcarrierIndices.^2));
    exponentialPhase = exp(1j * phaseTerm).';
    
    % 狄利克雷核
    thetaTerm = subcarrierIndices - q + dopplerShift + 2*numSubcarriers*chirpParam1*delayTap;
    
    numeratorTerm = exp(-1j * 2 * pi * thetaTerm) - 1;
    denominatorTerm = exp(-1j * 2 * pi/numSubcarriers * thetaTerm) - 1;
    
    dirichletKernel = zeros(numSubcarriers, 1);
    tol = 1e-6;
    
    % 处理分母为0的奇点
    singularityIndices = abs(denominatorTerm) < tol;
    dirichletKernel(~singularityIndices) = numeratorTerm(~singularityIndices) ./ denominatorTerm(~singularityIndices);
    dirichletKernel(singularityIndices) = numSubcarriers; 
    
    responseVector = (1/numSubcarriers) * exponentialPhase .* dirichletKernel;
    
end
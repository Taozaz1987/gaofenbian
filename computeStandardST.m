function [stTFR, tfrTime, tfrFreq] = computeStandardST(signal, t)
    % 计算标准S变换
    N = length(signal);
    dt = t(2) - t(1);
    
    % 频率范围
    fmin = 1; % Hz
    fmax = 100; % Hz
    fstep = 1; % Hz
    freqs = fmin:fstep:fmax;
    nFreqs = length(freqs);
    
    % 初始化时频矩阵
    stTFR = zeros(nFreqs, N);
    sigma=zeros(size(freqs));
    for fIdx = 1:nFreqs
        f = freqs(fIdx);
        sigma(fIdx) = 1 / abs(f);
        for tau = 1:N
            % 标准S变换窗口（高斯窗口）
           
            window = 1/(sigma(fIdx)*sqrt(2*pi)) * exp(-0.5 * ((t - t(tau))/sigma(fIdx)).^2);
            
            % 计算时频系数
            integrand = signal .* window .* exp(-1i * 2 * pi * f * t);
            stTFR(fIdx, tau) = trapz(t, integrand);
        end
    end
    
    tfrTime = t;
    tfrFreq = freqs;
end
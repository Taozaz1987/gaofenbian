%% 生成单个Ricker子波的函数
function w = ricker_wavelet(f0, t0, amplitude, t)
    % 生成零相位雷克子波
    % f0: 主频, t0: 峰值时间, amplitude: 振幅, t: 时间向量
    
    % 计算子波
    tau = pi * f0 * (t - t0);
    w = amplitude * (1 - 2 * tau.^2) .* exp(-tau.^2);
end

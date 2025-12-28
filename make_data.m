%% 1. 参数设置
clear; clc; close all;

% 基本参数
fs = 500;           % 采样率 (Hz)
dt = 1/fs;          % 采样间隔 (s)
T = 2.0;            % 记录长度 (s)
t = 0:dt:T-dt;      % 时间轴
N = length(t);      % 采样点数

% Q衰减参数
Q = 50;             % 品质因子（衰减强度，Q越小衰减越强）
f0 = 30;            % 子波主频 (Hz)

%% 2. 生成反射系数序列
% 方法1：随机反射系数（白噪声）
% r = randn(1, N) * 0.1;  % 高斯白噪声

% 方法2：几个孤立的反射界面（更真实）
r = zeros(1, N);
% 设置反射界面位置（秒）
reflector_times = [0.3, 0.6, 0.9, 1.2, 1.5];
% 设置反射系数大小（正负表示阻抗增加或减少）
reflector_amplitudes = [0.8, -0.6, 0.4, -0.3, 0.5];

for i = 1:length(reflector_times)
    idx = round(reflector_times(i) * fs) + 1;
    if idx <= N
        r(idx) = reflector_amplitudes(i);
    end
end

% 添加一些随机噪声使更真实
r = r + randn(1, N) * 0.02;

figure('Position', [100, 100, 800, 400]);
subplot(2,1,1);
plot(t, r, 'b-', 'LineWidth', 1.5);
xlabel('时间 (s)'); ylabel('振幅');
title('反射系数序列');
grid on;

%% 3. 生成地震子波
% Ricker子波（零相位）
t_wavelet = -0.1:dt:0.1;  % 子波时间轴
nt_wavelet = length(t_wavelet);
t_wavelet_center = floor(nt_wavelet/2) + 1;

% Ricker子波公式
ricker = (1 - 2*(pi*f0*t_wavelet).^2) .* exp(-(pi*f0*t_wavelet).^2);

% 确保子波能量归一化
ricker = ricker / max(abs(ricker));

subplot(2,1,2);
plot(t_wavelet, ricker, 'r-', 'LineWidth', 1.5);
xlabel('时间 (s)'); ylabel('振幅');
title(sprintf('Ricker子波 (主频 %.1f Hz)', f0));
grid on;

%% 4. 生成未衰减的地震记录（稳态卷积）
% 使用卷积生成地震记录
s_ideal = conv(r, ricker, 'same');  % 'same'保持相同长度

% 归一化
s_ideal = s_ideal / max(abs(s_ideal));

figure('Position', [100, 100, 1200, 400]);
subplot(1,3,1);
plot(t, s_ideal, 'b-', 'LineWidth', 1.5);
xlabel('时间 (s)'); ylabel('振幅');
title('未衰减的地震记录');
grid on;
xlim([0, T]);

%% 5. 应用Q衰减（时变滤波）
% 方法：将地震记录分成多个时窗，对每个时窗应用不同的衰减

% 定义衰减函数（常数Q模型）
% 衰减幅度：A(f, t) = exp(-π * f * t / Q)
% 同时考虑最小相位延迟

% 频率轴
f = (0:N/2) * fs / N;  % 正频率部分
f = f(2:end);          % 去掉直流分量
Nf = length(f);

% 初始化衰减后的信号
s_attenuated = zeros(1, N);

% 将信号分成多个重叠的时窗
window_length = round(0.2 * fs);  % 窗长度（200 ms）
overlap = round(0.5 * window_length);  % 50%重叠
step = window_length - overlap;

% 汉宁窗
hann_win = hanning(window_length)';

num_windows = floor((N - window_length) / step) + 1;

for win_idx = 1:num_windows
    % 当前时窗的起始和结束索引
    start_idx = (win_idx-1) * step + 1;
    end_idx = start_idx + window_length - 1;
    
    if end_idx > N
        end_idx = N;
    end
    
    % 提取当前时窗的信号
    win_signal = s_ideal(start_idx:end_idx);
    
    % 应用窗函数
    win_signal_windowed = win_signal .* hann_win(1:length(win_signal));
    
    % 计算当前时窗中心时间（旅行时）
    t_center = t(round((start_idx+end_idx)/2));
    
    % FFT
    win_fft = fft(win_signal_windowed);
    win_fft = win_fft(1:length(win_signal_windowed));  % 取前半部分
    
    % 计算衰减因子
    % 幅度衰减：exp(-π * f * t_center / Q)
    % 相位延迟：exp(-1i * π * f * t_center / Q) （最小相位近似）
    n_fft = length(win_fft);
    f_win = (0:n_fft/2) * fs / n_fft;
    f_win = f_win(2:end);  % 去掉直流
    
    % 创建衰减滤波器
    n_freq = length(f_win);
    attenuation_filter = zeros(1, n_fft);
    
    % 正频率部分
    for i = 1:n_freq
        % 幅度衰减
        amp_atten = exp(-pi * f_win(i) * t_center / Q);
        % 相位延迟（最小相位）
        phase_shift = exp(-1i * pi * f_win(i) * t_center / Q);
        % 综合衰减
        attenuation_filter(i+1) = amp_atten * phase_shift;
    end
    
    % 负频率部分（共轭对称）
    attenuation_filter(end-n_freq+1:end) = conj(attenuation_filter(n_freq+1:-1:2));
    
    % 应用衰减滤波器
    win_fft_attenuated = win_fft .* attenuation_filter(1:length(win_fft));
    
    % 逆FFT
    win_signal_attenuated = real(ifft(win_fft_attenuated));
    
    % 去窗并叠加到输出信号
    s_attenuated(start_idx:end_idx) = s_attenuated(start_idx:end_idx) + ...
        win_signal_attenuated .* hann_win(1:length(win_signal));
end

% 归一化
s_attenuated = s_attenuated / max(abs(s_attenuated));

subplot(1,3,2);
plot(t, s_attenuated, 'r-', 'LineWidth', 1.5);
xlabel('时间 (s)'); ylabel('振幅');
title(sprintf('Q衰减后的地震记录 (Q=%d)', Q));
grid on;
xlim([0, T]);

%% 6. 更精确的Q衰减方法（逐点衰减）
% 这种方法更精确但计算量更大
% 使用非稳态卷积模型

s_attenuated2 = zeros(1, N);

% 对每个反射系数应用时变衰减
for i = 1:N
    if abs(r(i)) > 0.01  % 只处理显著的反射系数
        % 该反射系数的到达时间
        t_arrival = t(i);
        
        % 计算衰减后的子波
        % 子波幅度衰减：exp(-π * f0 * t_arrival / Q)
        amp_atten = exp(-pi * f0 * t_arrival / Q);
        
        % 生成衰减后的子波（幅度减小，主频降低）
        f0_atten = f0 * amp_atten;  % 近似：频率也衰减
        ricker_atten = (1 - 2*(pi*f0_atten*t_wavelet).^2) .* exp(-(pi*f0_atten*t_wavelet).^2);
        ricker_atten = ricker_atten * amp_atten;
        
        % 将衰减子波放置在正确位置
        start_idx = i - floor(length(t_wavelet)/2);
        end_idx = start_idx + length(t_wavelet) - 1;
        
        % 确保索引在范围内
        if start_idx < 1
            ricker_part = ricker_atten(2-start_idx:end);
            start_idx = 1;
        else
            ricker_part = ricker_atten;
        end
        
        if end_idx > N
            ricker_part = ricker_part(1:end-(end_idx-N));
            end_idx = N;
        end
        
        % 添加到地震记录
        valid_length = end_idx - start_idx + 1;
        s_attenuated2(start_idx:end_idx) = s_attenuated2(start_idx:end_idx) + ...
            r(i) * ricker_part(1:valid_length);
    end
end

% 归一化
s_attenuated2 = s_attenuated2 / max(abs(s_attenuated2));

subplot(1,3,3);
plot(t, s_attenuated2, 'g-', 'LineWidth', 1.5);
xlabel('时间 (s)'); ylabel('振幅');
title(sprintf('逐点Q衰减 (Q=%d)', Q));
grid on;
xlim([0, T]);

%% 7. 频谱分析
figure('Position', [100, 100, 1000, 600]);

% 计算频谱
[P_ideal, f_axis] = pwelch(s_ideal, hamming(256), 128, 256, fs);
[P_atten, ~] = pwelch(s_attenuated, hamming(256), 128, 256, fs);
[P_atten2, ~] = pwelch(s_attenuated2, hamming(256), 128, 256, fs);

% 转换为dB
P_ideal_db = 10*log10(P_ideal);
P_atten_db = 10*log10(P_atten);
P_atten2_db = 10*log10(P_atten2);

subplot(2,2,1);
plot(t, s_ideal, 'b-', 'LineWidth', 1.5); hold on;
plot(t, s_attenuated, 'r-', 'LineWidth', 1.5);
xlabel('时间 (s)'); ylabel('振幅');
title('时域信号对比');
legend('未衰减', '衰减后', 'Location', 'northeast');
grid on;
xlim([0, T]);

subplot(2,2,2);
plot(f_axis, P_ideal_db, 'b-', 'LineWidth', 1.5); hold on;
plot(f_axis, P_atten_db, 'r-', 'LineWidth', 1.5);
xlabel('频率 (Hz)'); ylabel('功率谱密度 (dB)');
title('频谱对比');
legend('未衰减', '衰减后', 'Location', 'northeast');
grid on;
xlim([0, 100]);

subplot(2,2,3);
% 时频谱 - 未衰减
[s_ideal_tf, f_tf, t_tf] = spectrogram(s_ideal, 128, 120, 128, fs);
imagesc(t_tf, f_tf, 20*log10(abs(s_ideal_tf)));
axis xy; colorbar;
xlabel('时间 (s)'); ylabel('频率 (Hz)');
title('未衰减信号的时频谱');
ylim([0, 100]);

subplot(2,2,4);
% 时频谱 - 衰减后
[s_atten_tf, f_tf, t_tf] = spectrogram(s_attenuated, 128, 120, 128, fs);
imagesc(t_tf, f_tf, 20*log10(abs(s_atten_tf)));
axis xy; colorbar;
xlabel('时间 (s)'); ylabel('频率 (Hz)');
title('衰减后信号的时频谱');
ylim([0, 100]);

colormap(jet);

%% 8. 保存数据
% 可以选择保存生成的数据用于后续处理
data.t = t;
data.fs = fs;
data.r = r;                    % 反射系数
data.wavelet = ricker;         % 子波
data.s_ideal = s_ideal;        % 未衰减信号
data.s_attenuated = s_attenuated;  % 衰减信号
data.Q = Q;
data.f0 = f0;

save('synthetic_seismic_data.mat', 'data');

fprintf('地震记录生成完成！\n');
fprintf('参数：\n');
fprintf('  采样率：%.0f Hz\n', fs);
fprintf('  记录长度：%.1f s\n', T);
fprintf('  子波主频：%.0f Hz\n', f0);
fprintf('  品质因子 Q：%d\n', Q);
fprintf('  数据已保存到 synthetic_seismic_data.mat\n');

%% 9. 可选：添加随机噪声
% 如果需要更真实的数据，可以添加随机噪声
snr_db = 20;  % 信噪比 (dB)
noise_power = var(s_attenuated) / (10^(snr_db/10));
noise = sqrt(noise_power) * randn(1, N);
s_noisy = s_attenuated + noise;

figure('Position', [100, 100, 800, 400]);
subplot(1,2,1);
plot(t, s_attenuated, 'b-', 'LineWidth', 1.5); hold on;
plot(t, s_noisy, 'r-', 'LineWidth', 1);
xlabel('时间 (s)'); ylabel('振幅');
title(sprintf('添加噪声 (SNR = %.0f dB)', snr_db));
legend('无噪声', '有噪声');
grid on;
xlim([0, T]);

subplot(1,2,2);
[P_clean, f_axis] = pwelch(s_attenuated, hamming(256), 128, 256, fs);
[P_noisy, ~] = pwelch(s_noisy, hamming(256), 128, 256, fs);
plot(f_axis, 10*log10(P_clean), 'b-', 'LineWidth', 1.5); hold on;
plot(f_axis, 10*log10(P_noisy), 'r-', 'LineWidth', 1);
xlabel('频率 (Hz)'); ylabel('功率谱密度 (dB)');
title('频谱对比');
legend('无噪声', '有噪声');
grid on;
xlim([0, 100]);
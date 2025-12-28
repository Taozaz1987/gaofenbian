clear;clc;
%% 合成地震信号 - 基于ALFMT的动态反褶积方法验证

%% 1. 参数设置
fs = 1000;               % 采样频率 500 Hz
dt = 1/fs;              % 采样间隔 0.002 s
T = 1;                  % 记录长度 1 秒
t = 0:dt:T-dt;          % 时间轴
nt = length(t);         % 时间点数

% 质量因子 Q
Q = 80;

% 主频
f0 = 30;                % 30 Hz

%% 2. 生成最小相位子波
% 方法：先生成零相位雷克子波，然后转换为最小相位
% 生成零相位雷克子波
wavelet_length = 0.2;   % 子波长度 200ms
nw = round(wavelet_length * fs);
tw = (-nw/2:nw/2) * dt;

% 雷克子波公式
ricker_zero = (1 - 2*(pi*f0*tw).^2) .* exp(-(pi*f0*tw).^2);

% 将零相位子波转换为最小相位
% 方法：计算子波的希尔伯特变换，取其实部
ricker_min = real(hilbert(ricker_zero));

% 归一化
ricker_min = ricker_min / max(abs(ricker_min));

% 调整子波位置，使其在时间轴中心
center_idx = floor(length(t)/2);
start_idx = center_idx - floor(length(ricker_min)/2);
end_idx = start_idx + length(ricker_min) - 1;

% 确保索引在范围内
if start_idx < 1
    start_idx = 1;
end
if end_idx > nt
    end_idx = nt;
end

% 创建完整子波向量
w = zeros(1, nt);
w(start_idx:end_idx) = ricker_min(1:(end_idx-start_idx+1));

% 绘制子波
figure('Position', [100, 100, 800, 400]);
subplot(1,2,1);
plot(tw*1000, ricker_zero, 'b-', 'LineWidth', 1.5);
hold on;
plot(tw*1000, ricker_min, 'r-', 'LineWidth', 1.5);
title('雷克子波');
xlabel('时间 (ms)'); ylabel('振幅');
legend('零相位', '最小相位');
grid on;

subplot(1,2,2);
plot(t, w, 'k-', 'LineWidth', 1.5);
title('最小相位子波（时间域）');
xlabel('时间 (s)'); ylabel('振幅');
grid on;

%% 3. 生成反射系数序列（高斯白噪声）
% 设置随机种子以确保可重复性
rng(42);

% 生成高斯白噪声
reflection_coef = randn(1, nt);

% 对反射系数进行平滑处理，模拟地层连续性
window_size = 5;
smooth_filter = ones(1, window_size)/window_size;
reflection_coef = conv(reflection_coef, smooth_filter, 'same');

% 归一化反射系数
reflection_coef = reflection_coef / max(abs(reflection_coef)) * 0.5;

% 绘制反射系数
figure('Position', [100, 100, 800, 300]);
plot(t, reflection_coef, 'b-', 'LineWidth', 1);
title('反射系数序列（高斯白噪声）');
xlabel('时间 (s)'); ylabel('反射系数');
grid on;

%% 4. 生成衰减函数（基于Q模型）
% 根据常Q衰减理论：A(f,τ) = exp(-πfτ/Q)
% 其中：f是频率，τ是旅行时

% 计算频率轴
nfft = 2^nextpow2(nt);
f = (0:nfft-1) * (fs/nfft);
f_pos = f(1:nfft/2+1);  % 正频率部分

% 初始化衰减函数矩阵
attenuation_matrix = zeros(length(f_pos), nt);

% 计算每个时间点τ和每个频率f的衰减
for tau_idx = 1:nt
    tau = t(tau_idx);  % 旅行时
    
    % 计算振幅衰减
    amp_attenuation = exp(-pi * f_pos * tau / Q);
    
    % 计算相位畸变（希尔伯特变换关系）
    % 根据论文公式(20): a_Q(τ,f) = exp(-πfτ/Q + iH(πfτ/Q))
    % 其中H表示希尔伯特变换
    phase_term = pi * f_pos * tau / Q;
    
    % 计算希尔伯特变换（使用解析信号方法）
    analytic_signal = hilbert(phase_term);
    phase_attenuation = imag(analytic_signal);
    
    % 构建复数衰减函数
    attenuation_matrix(:, tau_idx) = amp_attenuation .* exp(1i * phase_attenuation);
end

%% 5. 生成非平稳地震信号
% 根据论文公式(21): s(f) = w(f) ∫ a_Q(τ,f) r(τ) exp(i2πfτ) dτ

% 计算子波频谱
W = fft(w, nfft);
W_pos = W(1:nfft/2+1);

% 初始化频谱
S_nonstationary = zeros(1, nfft);

% 对每个频率计算积分
for f_idx = 1:length(f_pos)
    f_current = f_pos(f_idx);
    
    % 计算积分项: ∫ a_Q(τ,f) r(τ) exp(i2πfτ) dτ
    integral_term = 0;
    for tau_idx = 1:nt
        tau = t(tau_idx);
        integral_term = integral_term + ...
                       attenuation_matrix(f_idx, tau_idx) * ...
                       reflection_coef(tau_idx) * ...
                       exp(1i * 2 * pi * f_current * tau) * dt;
    end
    
    % 计算该频率的频谱值
    S_nonstationary(f_idx) = W_pos(f_idx) * integral_term;
    
    % 设置负频率（共轭对称）
    if f_idx > 1 && f_idx < length(f_pos)
        S_nonstationary(nfft - f_idx + 2) = conj(S_nonstationary(f_idx));
    end
end

% 逆傅里叶变换得到时间域信号
s_nonstationary = real(ifft(S_nonstationary, nfft));
s_nonstationary = s_nonstationary(1:nt);  % 截取有效长度

% 归一化
s_nonstationary = s_nonstationary / max(abs(s_nonstationary));

%% 6. 作为对比，生成平稳地震信号（无衰减）
% 直接卷积：s(t) = w(t) * r(t)
s_stationary = conv(reflection_coef, w, 'same');
s_stationary = s_stationary / max(abs(s_stationary));

%% 7. 可视化结果
figure('Position', [100, 100, 1200, 800]);

% 子图1：子波和反射系数
subplot(3, 3, [1, 4]);
plot(t, w, 'b-', 'LineWidth', 1.5);
hold on;
plot(t, reflection_coef*0.2, 'r-', 'LineWidth', 1); % 缩放显示
title('最小相位子波和反射系数');
xlabel('时间 (s)'); ylabel('振幅');
legend('子波', '反射系数（缩放）', 'Location', 'northwest');
grid on;
xlim([0, 1]);

% 子图2：平稳地震信号
subplot(3, 3, 2);
plot(t, s_stationary, 'b-', 'LineWidth', 1.5);
title('平稳地震信号（无衰减）');
xlabel('时间 (s)'); ylabel('振幅');
grid on;
xlim([0, 1]);

% 子图3：非平稳地震信号
subplot(3, 3, 3);
plot(t, s_nonstationary, 'r-', 'LineWidth', 1.5);
title(sprintf('非平稳地震信号（Q=%d）', Q));
xlabel('时间 (s)'); ylabel('振幅');
grid on;
xlim([0, 1]);

% 子图4：对比图
subplot(3, 3, [5, 6]);
plot(t, s_stationary, 'b-', 'LineWidth', 1.5, 'DisplayName', '平稳信号');
hold on;
plot(t, s_nonstationary, 'r-', 'LineWidth', 1.5, 'DisplayName', sprintf('非平稳信号（Q=%d）', Q));
title('平稳与非平稳信号对比');
xlabel('时间 (s)'); ylabel('振幅');
legend('Location', 'best');
grid on;
xlim([0, 1]);

% 子图5：衰减函数示例（某个时间点的频率响应）
subplot(3, 3, 7);
tau_example = 0.5; % 0.5秒处的衰减
tau_idx = find(t >= tau_example, 1);
plot(f_pos, abs(attenuation_matrix(:, tau_idx)), 'k-', 'LineWidth', 1.5);
title(sprintf('衰减函数 (τ=%.1fs, Q=%d)', tau_example, Q));
xlabel('频率 (Hz)'); ylabel('|A(f,τ)|');
grid on;
xlim([0, 100]);

% 子图6：频谱对比
subplot(3, 3, 8);
[P_stationary, f_stationary] = pwelch(s_stationary, [], [], [], fs);
[P_nonstationary, f_nonstationary] = pwelch(s_nonstationary, [], [], [], fs);
plot(f_stationary, 10*log10(P_stationary), 'b-', 'LineWidth', 1.5, 'DisplayName', '平稳信号');
hold on;
plot(f_nonstationary, 10*log10(P_nonstationary), 'r-', 'LineWidth', 1.5, 'DisplayName', '非平稳信号');
title('频谱对比');
xlabel('频率 (Hz)'); ylabel('功率谱密度 (dB)');
legend('Location', 'best');
grid on;
xlim([0, 100]);

% 子图7：时频分析（STFT对比）
subplot(3, 3, 9);
% 平稳信号的STFT
[S_stationary, F_st, T_st] = spectrogram(s_stationary, 64, 60, 128, fs);
% 非平稳信号的STFT
[S_nonstationary_tf, F_nst, T_nst] = spectrogram(s_nonstationary, 64, 60, 128, fs);

% 计算平均频率
mean_freq_st = mean(abs(S_stationary), 1);
mean_freq_nst = mean(abs(S_nonstationary_tf), 1);

plot(T_st, mean_freq_st, 'b-', 'LineWidth', 1.5, 'DisplayName', '平稳信号');
hold on;
plot(T_nst, mean_freq_nst, 'r-', 'LineWidth', 1.5, 'DisplayName', '非平稳信号');
title('平均瞬时频率');
xlabel('时间 (s)'); ylabel('平均频率 (Hz)');
legend('Location', 'best');
grid on;
xlim([0, 1]);

%% 8. 输出统计信息
fprintf('========== 合成地震信号统计信息 ==========\n');
fprintf('采样频率: %d Hz\n', fs);
fprintf('记录长度: %.1f s\n', T);
fprintf('质量因子 Q: %d\n', Q);
fprintf('子波主频: %d Hz\n', f0);
fprintf('信号长度: %d 点\n', nt);
fprintf('平稳信号最大振幅: %.4f\n', max(abs(s_stationary)));
fprintf('非平稳信号最大振幅: %.4f\n', max(abs(s_nonstationary)));
fprintf('能量衰减比: %.2f%%\n', ...
    100 * sum(s_nonstationary.^2) / sum(s_stationary.^2));

% 计算主频偏移
[~, idx_st] = max(abs(fft(s_stationary, nfft)));
[~, idx_nst] = max(abs(fft(s_nonstationary, nfft)));
fprintf('平稳信号主频: %.1f Hz\n', f(idx_st));
fprintf('非平稳信号主频: %.1f Hz\n', f(idx_nst));

%% 9. 保存数据
% 可以选择保存数据用于后续处理
%save('synthetic_seismic_data.mat', ...
%     't', 's_stationary', 's_nonstationary', ...
%     'reflection_coef', 'w', 'Q', 'f0', 'fs');
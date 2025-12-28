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
s=s_nonstationary;
%% ALFMT
k=4;%假设k=
f=linspace(0,100-100/1000,1000);
f_optloc=optloc(s,t);%计算局部频率
ALFMT=zeros(length(t),length(f));
dt=t(2)-t(1);
df=f(2)-f(1);
sigma=zeros(length(t),length(f));
for tau_idx = 1:length(t)
    tau = t(tau_idx);
    f_loc_val = f_optloc(tau_idx);  % 获取该时刻的局部频率
    
    % 向量化计算sigma（针对所有频率点）
    sigma_inv = (f_loc_val + abs(f_loc_val - f)) / k;  % 向量
    sigma(tau_idx,:) = 1 ./ sigma_inv;  % 实际的sigma值
    
    % 预计算高斯窗和频率项
    for f_idx = 1:length(f)
        % 高斯窗（向量化）
        window = exp(-(t - tau).^2 .* (sigma_inv(f_idx)^2)/2);
        
        % 完整被积函数
        integrand = s .* window .* exp(-1i * 2 * pi * f(f_idx) * t);
        
        % 使用trapz积分
        ALFMT(tau_idx, f_idx) = sigma_inv(f_idx)/sqrt(2*pi) * trapz(t, integrand);
    end
end
magnitude_ALFMT=abs(ALFMT);
ALFMT_s_hs=hs(ALFMT,t,f);
phi=calculate_phi(ALFMT_s_hs,t);
V_ALFMT_r=calculate_Vr(ALFMT,ALFMT_s_hs,phi);
r_reconstruct=de_V_ALFMT(V_ALFMT_r,t,f);
r_real=real(r_reconstruct);
s_recon=conv(r_real,ricker_min,'same');

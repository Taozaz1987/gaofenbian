clear;clc;
%% 信号参数设置
fs = 1000;          % 采样频率 (Hz)
T = 1;              % 信号时长 (秒)
t = 0:1/fs:T-1/fs;  % 时间向量
N = length(t);      % 采样点数

%% 1. 固定频率信号分量
A1 = 1.0; f1 = 10;   % 10Hz 信号
A2 = 0.8; f2 = 30;   % 30Hz 信号  
A3 = 0.6; f3 = 60;   % 60Hz 信号
A4 = 0.4; f4 = 80;   % 80Hz 信号

% 生成正弦信号
signal_f1 = A1 * sin(2*pi*f1*t);
signal_f2 = A2 * sin(2*pi*f2*t);
signal_f3 = A3 * sin(2*pi*f3*t);
signal_f4 = A4 * sin(2*pi*f4*t);

%% 2. 调频信号 (FM Signal)
% 调频信号参数
A_fm = 1.2;          % FM信号幅值
%f_carrier = 40;      % 载波频率 (Hz)
%f_modulation = 2;    % 调制频率 (Hz)
%modulation_index = 5; % 调制指数
f_ch=20+50*t;
% 方法1: 使用sin(2*pi*积分(瞬时频率))生成FM信号
% 瞬时频率: f_instant = f_carrier + modulation_index*f_modulation*cos(2*pi*f_modulation*t)
%phase_integral = cumsum(2*pi*(f_carrier + modulation_index*f_modulation*cos(2*pi*f_modulation*t)))/fs;
%signal_fm = A_fm * sin(phase_integral);

% 方法2: 简化的FM表达式 (对于正弦调制)
signal_fm = A_fm * sin(2*pi.*f_ch.*t);

%% 3. 组合所有信号
signal_total = signal_f1 + signal_f2 + signal_f3 + signal_f4 + signal_fm;

%% 4. 可选：添加高斯白噪声
SNR = 20;  % 信噪比 (dB)
signal_noisy = awgn(signal_total, SNR, 'measured');

%% 绘制时域信号
figure('Position', [100, 100, 1200, 800])

% 子图1: 各分量信号
subplot(4,2,1)
plot(t, signal_f1, 'b')
title('10Hz 信号')
xlabel('时间 (s)'); ylabel('幅值')
xlim([0, 0.5]); grid on

subplot(4,2,3)
plot(t, signal_f2, 'r')
title('30Hz 信号')
xlabel('时间 (s)'); ylabel('幅值')
xlim([0, 0.2]); grid on

subplot(4,2,5)
plot(t, signal_f3, 'g')
title('60Hz 信号')
xlabel('时间 (s)'); ylabel('幅值')
xlim([0, 0.1]); grid on

subplot(4,2,7)
plot(t, signal_f4, 'm')
title('80Hz 信号')
xlabel('时间 (s)'); ylabel('幅值')
xlim([0, 0.1]); grid on

% 子图2: FM信号
subplot(4,2,2)
plot(t, signal_fm, 'Color', [0.8, 0.4, 0])
title('调频信号 (40Hz载波，2Hz调制)')
xlabel('时间 (s)'); ylabel('幅值')
xlim([0, 1]); grid on

% 子图3: 组合信号
subplot(4,2,[4,6])
plot(t, signal_total, 'k', 'LineWidth', 1.2)
title('合成信号 (10+30+60+80Hz + FM)')
xlabel('时间 (s)'); ylabel('幅值')
xlim([0, 1]); grid on

% 子图4: 含噪声的信号
subplot(4,2,8)
plot(t, signal_noisy, 'Color', [0.3, 0.3, 0.3])
title(sprintf('添加噪声的信号 (SNR=%ddB)', SNR))
xlabel('时间 (s)'); ylabel('幅值')
xlim([0, 0.2]); grid on

sgtitle('多频率信号与调频信号合成')
s=signal_total;
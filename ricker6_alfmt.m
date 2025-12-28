clear;clc;
%% 参数设置
fs = 1000;          % 采样频率 1000 Hz
dt = 1/fs;          % 采样间隔 0.001 s
T = 1;              % 总时长 1 秒
t = 0:dt:T-dt;      % 时间轴

%% 定义6个子波的参数
% 每个子波的 [主频(Hz), 峰值时间(s), 振幅]
wavelet_params = [
    30, 0.10, 1.0;   % 20Hz 在 0.1s
    10, 0.20, 0.8;   % 10Hz 在 0.2s
    40, 0.35, -1.2;   % 20Hz 在 0.35s
    50, 0.50, 0.9;   % 30Hz 在 0.5s
    20, 0.65, -1.1;   % 20Hz 在 0.65s
    70, 0.85, 1.0;   % 30Hz 在 0.85s
];

%% 生成合成地震记录
seismic_trace = zeros(size(t));


% 绘制每个子波和叠加过程
for i = 1:size(wavelet_params, 1)
    f0 = wavelet_params(i, 1);
    t0 = wavelet_params(i, 2);
    amp = wavelet_params(i, 3);
    
    % 生成当前子波
    w = ricker_wavelet(f0, t0, amp, t);
    
    % 叠加到总地震记录
    seismic_trace = seismic_trace + w;
end
s=seismic_trace;
%% ALFMT
k=1.2;%假设k=
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
figure(1);
p1=plot(t,s,'b-');
figure(2);
imagesc(t, f, magnitude_ALFMT');
axis xy;
colorbar;
figure(3);
[stTFR, tfrTime, tfrFreq] = computeStandardST(s, t);
imagesc(tfrTime, tfrFreq, abs(stTFR));
axis xy;
colorbar;

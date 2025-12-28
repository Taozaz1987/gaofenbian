%% 测试信号
t = linspace(0,1,1000);
s = 2.*cos(2.*pi.*(80.*t)) + ...
    (1+0.5.*cos(2.*t)).*exp(-t/10).*cos(10.*pi.*(8.*t+6.*t.^2) + 0.3.*cos(t)) + ...
    (2+0.2.*cos(t)).*sin(10*pi.*(5.*t+0.3.*cos(6.*t)));
figure(1); plot(t, s, 'b', 'LineWidth', 1.5); hold on;

%% 参数设置
k = 1;  % ALFMT参数
dt = t(2)-t(1);
f = linspace(0,100,1000); % 频率轴
df = f(2)-f(1);

%% 计算局部频率
f_optloc = optloc(s, t);

%% ALFMT正变换 (修正高斯窗)
ALFMT = zeros(length(t), length(f));
for tau = 1:length(t)
    for f_index = 1:length(f)
        % 修正1: 正确计算sigma
        sigma = (f_optloc(tau) + abs(f_optloc(tau) - f(f_index))) / k;
        % 修正2: 高斯窗指数正确形式
        gauss_win = exp(-(t - t(tau)).^2 / (2 * sigma^2)); % 向量化计算
        % 修正3: 归一化因子 (保持物理意义)
        zuo = 1 / (sigma * sqrt(2*pi)); % k=1时简化
        % 向量化积分 (避免内层循环)
        integrand = s .* gauss_win .* exp(-1i*2*pi*f(f_index)*t);
        you = trapz(t, integrand); % 使用梯形积分更准确
        ALFMT(tau, f_index) = zuo * you;
    end
end

%% ALFMT逆变换 (添加窗函数补偿)
x_rec = zeros(size(t));
C = 0; % 归一化常数

for t_index = 1:length(t)
    for tau = 1:length(t)
        for f_index = 1:length(f)
            % 重新计算当前(tau,f)处的sigma
            sigma = (f_optloc(tau) + abs(f_optloc(tau) - f(f_index))) / k;
            % 修正4: 逆变换必须包含相同的高斯窗
            g_tau = exp(-(t(t_index) - t(tau))^2 / (2 * sigma^2));
            % 累加重建信号
            x_rec(t_index) = x_rec(t_index) + ...
                ALFMT(tau, f_index) * g_tau * exp(1i*2*pi*f(f_index)*t(t_index)) * df;
        end
        % 累加归一化常数 (框架界)
        C = C + exp(-(t(t_index) - t(tau))^2 / (2 * sigma^2))^2 * dt;
    end
    x_rec(t_index) = x_rec(t_index) * dt;
end

% 修正5: 全局归一化 (解决边界能量损失)
x_rec = real(x_rec) / max(C); 

%% 绘制结果
plot(t, x_rec, 'r--', 'LineWidth', 1.5);
legend('Original Signal', 'Reconstructed Signal');
title('ALFMT Reversibility Verification');
xlabel('Time (s)'); ylabel('Amplitude');
grid on;

%% 局部频率计算 (保持原函数，但优化数值稳定性)
function f_optloc = optloc(x, t)
    x_c = hilbert(x);
    h = imag(x_c);
    h_dao = gradient(h, t);
    x_dao = gradient(x, t);
    n = x .* h_dao - x_dao .* h;
    D = 2*pi*(x.^2 + h.^2);
    epsilon = max(abs(D)) * 1e-5; % 自适应正则化
    f_optloc = n ./ (D + epsilon); % 逐元素除法
end
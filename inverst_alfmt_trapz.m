clear;clc;
%% 测试信号
fs = 1000;
t = 0:1/fs:1-1/fs;
s=2.*cos(2.*pi.*(80.*t))+(1+0.5.*cos(2.*t)).*exp(-t/10).*cos(10.*pi.*(8.*t+6.*t.^2)...
    +0.3.*cos(t))+(2+0.2.*cos(t)).*sin(10*pi.*(5.*t+0.3.*cos(6.*t)));

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

%% inv_ALFMT
x=zeros(1,length(s));
for t_index=1:length(t )
    %f_loc = f_optloc(t_index);
    for f_index=1:length(f )
        fk = f(f_index);
        for tau=1:length(t )
            x(t_index)=ALFMT(tau,f_index)*exp(1i*2*pi*fk...
                *t(t_index))*dt+x(t_index);
        end        
    end
    x(t_index)=x(t_index)*df;
end
real_x=real(x);
test_x=2*real_x;%不知道为啥振幅对不上

figure(1);
p1=plot(t,s,'b-');
hold on;
p2=plot(t,test_x,'r--');
legend([p1, p2], {'原始信号', '重构信号'}, ...
       'Location', 'northeast', ...  % 右上角
       'FontSize', 12);
xlabel('时间');
ylabel('幅值');
title('信号对比图');
grid on;
figure(2);
imagesc(t, f, magnitude_ALFMT');
axis xy;
colorbar;
xlabel('Time (s)');
ylabel('Frequency (Hz)');
title(['ALFMT Time-Frequency Spectrum K = ', num2str(k)]);
figure(3);
plot(t,f_optloc);





       



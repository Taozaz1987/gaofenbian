clear;clc;
%% 测试信号
t=linspace(0,1,1000);
%s=2*sin(2*pi*20*t)+1.5*cos(sin(pi*5*t));
s=2.*cos(2.*pi.*(80.*t))+(1+0.5.*cos(2.*t)).*exp(-t/10).*cos(10.*pi.*(8.*t+6.*t.^2)...
    +0.3.*cos(t))+(2+0.2.*cos(t)).*sin(10*pi.*(5.*t+0.3.*cos(6.*t)));

%% ALFMT
k=1;%假设k=1
f=linspace(0,100,1000);
f_optloc=loc(s,t);%计算局部频率
ALFMT=zeros(length(t),length(f));
dt=t(2)-t(1);
df=f(2)-f(1);
for tau=1:length(t)
    for f_index=1:length(f)
        you=0;
        sigma=(f_optloc(tau)+abs(f_optloc(tau)-f(f_index)))/k;%sigma^-1
        zuo=sigma/sqrt(2*pi);
        for t_index=1:length(t)
            you=you+s(t_index).*exp(-(t(t_index)-t(tau))^2/2*sigma^2)*...
                exp(-1i*2*pi*f(f_index)*t(t_index))*dt;
        end
        ALFMT(tau,f_index)=zuo*you;
    end
end
magnitude_ALFMT=abs(ALFMT);

%% inv_ALFMT
x=zeros(1,length(s));
for t_index=1:length(t )
    for f_index=1:length(f )
        for tau=1:length(t )
            x(t_index)=ALFMT(tau,f_index)*exp(1i*2*pi*f(f_index)...
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
title('ALFMT Time-Frequency Spectrum');






       



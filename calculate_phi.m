function phi=calculate_phi(ALFMT_s_hs,t)
phi=zeros(size(ALFMT_s_hs));
for tau=1:length(t)
    log_hs_A=log(ALFMT_s_hs(tau,:)+eps);%ln|V_alfmt_s|hs
    analytic_signal = hilbert(log_hs_A);%构造解析信号
    phi(tau,:)=imag(analytic_signal);%取虚部——相位
end
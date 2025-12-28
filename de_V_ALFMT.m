function r_reconstruct=de_V_ALFMT(V_ALFMT_r,t,f)
dt=t(2)-t(1);
df=f(2)-f(1);
R_f=sum(V_ALFMT_r,1)*dt;
exp_matrix=exp(1i*2*pi*(f(:))*t);
r_reconstruct=2*real(R_f*exp_matrix)*df;



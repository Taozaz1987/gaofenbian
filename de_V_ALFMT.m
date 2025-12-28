function r_reconstruct=de_V_ALFMT(V_ALFMT_r,t,f)
dt=t(2)-t(1);
df=f(2)-f(1);
r_reconstruct=zeros(1,length(t));
for t_index=1:length(t )
    for f_index=1:length(f )
        fk = f(f_index);
        for tau=1:length(t )
            r_reconstruct(t_index)=V_ALFMT_r(tau,f_index)*exp(1i*2*pi*fk...
                *t(t_index))*dt+r_reconstruct(t_index);
        end        
    end
    r_reconstruct(t_index)=r_reconstruct(t_index)*df;
end



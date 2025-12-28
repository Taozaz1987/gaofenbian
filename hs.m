function ALFMT_s_hs=hs(V_ALFMT_s,t,f)
N=10;
A=abs(V_ALFMT_s);
ALFMT_s_hs=zeros(size(A));
c_max=length(t)*length(f);
delta_c=c_max/N;
I=zeros(length(t),length(f),N);
for n=1:N
    for tau=1:length(t)
       for f_index=1:length(f)
           if (tau*f_index>=(n-1)*delta_c)&&(tau*f_index<n*delta_c)
               I(tau,f_index,n)=1;
           else 
               I(tau,f_index,n)=0;
           end
       end
    end
end
a=zeros(1,N);
for n=1:N
    A_In=A.*I(:,:,n);
    up=sum(A_In,'all');
    down=sum(I(:,:,n),'all')+eps;
    a(n)=up/down;
end
for n=1:N
    ALFMT_s_hs=ALFMT_s_hs+a(n).*I(:,:,n);
end

    
    

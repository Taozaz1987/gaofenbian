%% inst瞬时频率
function f_optloc=inst(x,t)
f_optloc=zeros(length(t),1);
x_c=hilbert(x);
h=imag(x_c);
h_dao=gradient(h,t);
x_dao=gradient(x,t);
n=zeros(length(t),1);
D=zeros(length(t),length(t));
for i=1:length(t)
    n(i)=x(i)*h_dao(i)-x_dao(i)*h(i);
    D(i,i)=2*pi*((x(i)^2)+(h(i)^2))+eps;
end
I=eye(length(D));
epsilon=max(n)*(1e-5);
f_optloc=(D+epsilon*I)\n;
end

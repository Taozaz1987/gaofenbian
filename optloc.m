function f_optloc=optloc(x,t)
dt=t(2)-t(1);
N=length(x);
l_x=length(x);
l_window=221;
half_window=(l_window-1)/2;
xi=9;%待定
epsilon=1e-4;%待定
A=zeros(size(x));
x_c=hilbert(x);
h=imag(x_c);
I=eye(length(x));
for i=1:length(x)
    A(i)=sqrt(x(i)^2+h(i)^2);
end
A2_average=sum(A.^2)/l_x;
A2_aver_per=zeros(size(x));
for i=1:length(x)
    if i-half_window<1
        A2_aver_per(i)=mean(A(1:i*2-1).^2);
    elseif i+half_window>length(x)
        A2_aver_per(i)=mean(A((2*i-length(x)):length(x)).^2);
    else
        A2_aver_per(i)=mean(A(i-half_window:i+half_window).^2);
    end
end
P=zeros(N-1,N);
for i=1:N-1
    for j=1:N
        if j==i
            P(i,j)=sqrt(A2_average/A2_aver_per(i));
        elseif j==i+1
            P(i,j)=-sqrt(A2_average/A2_aver_per(j));
        else
            P(i,j)= 0;
        end
    end
end
lamda=sqrt(A2_aver_per+epsilon^2);
d=2*pi.*A(:).^2;
D=diag(d);
h_dao=gradient(h,dt);
x_dao=gradient(x,dt);
n=zeros(length(x),1);
for i=1:length(x)
    n(i)=x(i).*h_dao(i)-x_dao(i).*h(i);
end
S_inv=I+xi^2.*P.'*P;
S=inv(I+xi^2.*P.'*P);
mid=inv(diag(lamda.^2)+S_inv\(D-diag(lamda.^2)));
f_optloc=mid*S*n;
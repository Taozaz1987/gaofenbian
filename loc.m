%% optloc局部频率
function f_loc=loc(x,t)
xi=0.5;%平滑算子
lamda=0.01;
%epsilon=0.1;
N=length(t);
K=zeros(N-1,N);
for i=1:N-1
    for j=1:N
        if j==i
            K(i,j)=-1;
        elseif j==i+1
            K(i,j)=1;
        else
            K(i,j)=0;
        end
    end
end
I=eye(length(t));
A=I+xi^2*(K'*K);%
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
%Sn=(I+xi^2*K'*K)^(-1)n=result
result=A\n;
%f_loc=(lamda^2*I+S*(D-lamda^2*I))^(-1)*S*n
B=lamda^2*I+A\(D-lamda^2*I);
f_loc=B\result;
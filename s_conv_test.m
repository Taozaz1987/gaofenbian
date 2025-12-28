clear;clc;
f0=30;%hz
t0=0.08;%
fs=1000;
t=0:1/fs:1-1/fs;
r=zeros(1,length(t));
r(100)=0.8;
r(250)=0.6;
r(500)=-0.8;
r(600)=0.7;
r(800)=-0.5;
r(860)=-0.85;
ricker=mk_ricker(f0,t0,fs);
s=conv(r, ricker, 'same'); 


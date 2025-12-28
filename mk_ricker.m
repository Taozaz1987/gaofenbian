%% 零相位leike子波
function ricker=mk_ricker(f0,t0,fs,A)
dt=1/fs;
nw=round(t0*fs);%子波一侧采样点数
tw = (-nw/2 : nw/2) * dt;
ricker =A* (1 - 2*(pi*f0*tw).^2) .* exp(-(pi*f0*tw).^2);

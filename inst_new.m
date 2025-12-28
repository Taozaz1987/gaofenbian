%% inst瞬时频率
function f_optloc=inst_new(x,fs)
h = hilbert(x);
f_optloc = fs/(2*pi)*diff(unwrap(angle(h)));
f_optloc=[f_optloc,f_optloc(end)];
for i=1:length(f_optloc)
   if f_optloc(i)>=100
        f_optloc(i)=100;
    elseif f_optloc(i)<=0
        f_optloc(i)=0;
    end
end
 f_optloc = smooth(f_optloc, 5);
end
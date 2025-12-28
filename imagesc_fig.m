imagesc(tfrTime, tfrFreq, abs(stTFR));
axis xy;
xlabel('时间 (s)');
ylabel('频率 (Hz)');
title('标准S变换');
colorbar;
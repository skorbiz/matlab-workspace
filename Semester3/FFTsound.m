clear;
data = importdata('tone1.txt');
FTData = fft(data);


subplot(2,1,1); plot(abs(FTData)),grid;
title('Amplitude Spectrum');
xlabel('Frequency (Hz)');
ylabel('|data|');
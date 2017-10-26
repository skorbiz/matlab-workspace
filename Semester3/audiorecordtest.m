clear;
Fs = 16000;                  % Sampling frequency
T = 1/Fs;                    % Sample time
L = 10000;                    % Length of signal
t = (0:L-1)*T;               % Time vector
nBits = 16;                  % Number of bits used in sample  
nChanels = 1;                % Channels

recObj = audiorecorder(Fs,nBits,nChanels);
recordblocking(recObj, T*L);
theRec = getaudiodata(recObj);     

subplot(3,1,1),plot(Fs*t,theRec)
title('Input signal')
xlabel('time (milliseconds)')

Y = fft(theRec)/L;
f = Fs/2*linspace(0,1,L/2+1);
%magnitude = 20*log10(abs(Y(1:L/2+1)));
magnitude = 2*abs(Y(1:L/2+1));

subplot(3,1,2),plot(f,magnitude);grid;
title('Single-Sided Amplitude Spectrum')
xlabel('Frequency (Hz)')
ylabel('Magnitude')

PF1 = uint64(697*L/Fs+1);
PF2 = uint64(770*L/Fs+1);
PF3 = uint64(852*L/Fs+1);
PF4 = uint64(941*L/Fs+1);
PF5 = uint64(1209*L/Fs+1);
PF6 = uint64(1336*L/Fs+1);
PF7 = uint64(1477*L/Fs+1);
PF8 = uint64(1633*L/Fs+1);

peak = magnitude([PF1 PF2 PF3 PF4 PF5 PF6 PF7 PF8])

SNR1  = 20*log10((3*peak(4))/(peak(2)+peak(3)+peak(1)));
SNR2  = 20*log10((3*peak(5))/(peak(6)+peak(7)+peak(8)));

SNR = [SNR1 SNR2];

subplot(3,1,3),bar(SNR, 'grouped');grid;
title('SNR')
xlabel('Test nr.')
ylabel('SNR')



%Afspilning af lyd:
%Fs = 4000;      % Samples per second
%toneFreq = 3000;  % Tone frequency, in Hertz
%nSeconds = 4;   % Duration of the sound
%y = sin(linspace(0,nSeconds*toneFreq*2*pi,round(nSeconds*Fs)));
%sound(y,Fs);  % Play sound at sampling rate Fs
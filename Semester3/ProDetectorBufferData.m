%Data behandling til Detectoring af en tone med forskellige buffere:
%clear;

% 16 bit - 8000 Hz
% Buffer 256
% Buffer 512
% Buffer 1024
% Buffer 4096
% Naturligt større buffer giver bedre signal

%Extract data
d1 = DetectorTest0x2D16bit0x2D80000x2Dbuffer0x2D256';
d2 = DetectorTest0x2D16bit0x2D80000x2Dbuffer0x2D512';
d3 = DetectorTest0x2D16bit0x2D80000x2Dbuffer0x2D1024';
d4 = DetectorTest0x2D16bit0x2D80000x2Dbuffer0x2D4096';


%Data omkring samplinger
Fs = 8000;                    % Sampling frequency
T = 1/Fs;                    % Sample time
L1 = 256;                       % Length of signal
L2 = 512;                       % Length of signal
L3 = 1024;                      % Length of signal
L4 = 4096;                      % Length of signal

t1 = (0:L1-1)*T;               % Time vector
t2 = (0:L2-1)*T;               % Time vector
t3 = (0:L3-1)*T;               % Time vector
t4 = (0:L4-1)*T;               % Time vector

%Behandling af data
f1 = Fs/2*linspace(0,1,L1/2+1);
f2 = Fs/2*linspace(0,1,L2/2+1);
f3 = Fs/2*linspace(0,1,L3/2+1);
f4 = Fs/2*linspace(0,1,L4/2+1);

theRec = d1;
Y = fft(theRec)/L1;
mag1 = 2*abs(Y(1:L1/2+1));
magdB1 = 20*log10(abs(Y(1:L1/2+1)));

theRec = d2;
Y = fft(theRec)/L2;
mag2 = 2*abs(Y(1:L2/2+1));
magdB2 = 20*log10(abs(Y(1:L2/2+1)));

theRec = d3;
Y = fft(theRec)/L3;
mag3 = 2*abs(Y(1:L3/2+1));
magdB3 = 20*log10(abs(Y(1:L3/2+1)));

theRec = d4;
Y = fft(theRec)/L4;
mag4 = 2*abs(Y(1:L4/2+1));
magdB4 = 20*log10(abs(Y(1:L4/2+1)));


%Plot
figure(1)
subplot(2,2,1),plot(f1,mag1,f2,mag2,f3,mag3,f4,mag4);grid;axis([0 4100 0 10000]);
title('Single-Sided Amplitude Spectrum')
xlabel('Frequency (Hz)')
ylabel('Magnitude')
subplot(2,2,2),plot(f1,magdB1,f2,magdB2,f3,magdB3,f4,magdB4);grid;axis([0 4100 -20 100]);
title('Single-Sided Amplitude Spectrum')
xlabel('Frequency (Hz)')
ylabel('Magnitude')

%Udregning af SNR
L = L1;
PF1 = uint64(697*L/Fs+1);
PF2 = uint64(770*L/Fs+1);
PF3 = uint64(852*L/Fs+1);
PF4 = uint64(941*L/Fs+1);
PF5 = uint64(1209*L/Fs+1);
PF6 = uint64(1336*L/Fs+1);
PF7 = uint64(1477*L/Fs+1);
PF8 = uint64(1633*L/Fs+1);

peak = mag1([PF1 PF2 PF3 PF4 PF5 PF6 PF7 PF8]);
SNR1  = 20*log10((3*peak(1))/(peak(2)+peak(3)+peak(4)));
SNR2  = 20*log10((3*peak(5))/(peak(6)+peak(7)+peak(8)));

L = L2
PF1 = uint64(697*L/Fs+1);
PF2 = uint64(770*L/Fs+1);
PF3 = uint64(852*L/Fs+1);
PF4 = uint64(941*L/Fs+1);
PF5 = uint64(1209*L/Fs+1);
PF6 = uint64(1336*L/Fs+1);
PF7 = uint64(1477*L/Fs+1);
PF8 = uint64(1633*L/Fs+1);


peak = mag2([PF1 PF2 PF3 PF4 PF5 PF6 PF7 PF8]);
SNR3  = 20*log10((3*peak(1))/(peak(2)+peak(3)+peak(4)));
SNR4  = 20*log10((3*peak(5))/(peak(6)+peak(7)+peak(8)));

L = L3
PF1 = uint64(697*L/Fs+1);
PF2 = uint64(770*L/Fs+1);
PF3 = uint64(852*L/Fs+1);
PF4 = uint64(941*L/Fs+1);
PF5 = uint64(1209*L/Fs+1);
PF6 = uint64(1336*L/Fs+1);
PF7 = uint64(1477*L/Fs+1);
PF8 = uint64(1633*L/Fs+1);


peak = mag3([PF1 PF2 PF3 PF4 PF5 PF6 PF7 PF8]);
SNR5  = 20*log10((3*peak(1))/(peak(2)+peak(3)+peak(4)));
SNR6  = 20*log10((3*peak(5))/(peak(6)+peak(7)+peak(8)));

L = L4
PF1 = uint64(697*L/Fs+1);
PF2 = uint64(770*L/Fs+1);
PF3 = uint64(852*L/Fs+1);
PF4 = uint64(941*L/Fs+1);
PF5 = uint64(1209*L/Fs+1);
PF6 = uint64(1336*L/Fs+1);
PF7 = uint64(1477*L/Fs+1);
PF8 = uint64(1633*L/Fs+1);


peak = mag4([PF1 PF2 PF3 PF4 PF5 PF6 PF7 PF8]);
SNR7  = 20*log10((3*peak(1))/(peak(2)+peak(3)+peak(4)));
SNR8  = 20*log10((3*peak(5))/(peak(6)+peak(7)+peak(8)));

SNR = [SNR1 SNR2; SNR3 SNR4; SNR5 SNR6; SNR7 SNR8];

subplot(2,1,2),bar(SNR, 'grouped');grid;
title('SNR')
xlabel('Test')
ylabel('SNR')


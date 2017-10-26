%Data behandling til Detectoring af en tone tone:
%clear;

% 16 bit - 8000 Hz - 1000 Samples
% No window
% Hamming
% Blackman
% Blackman viduet er langt bedre end både intet vidue og hamming. Hamming
% er dog også langt bedre end at have intet vindue.


%Extract data
d1 = DetectorTest0x2D16bit0x2D8000';
d2 = DetectorTest0x2D16bit0x2D80000x2Dblackman0x2Dwindow';
d3 = DetectorTest0x2D16bit0x2D80000x2Dhamming0x2Dwindow';


%Data omkring samplinger
Fs = 8000;                    % Sampling frequency
T = 1/Fs;                    % Sample time
L = 1000;                      % Length of signal
t = (0:L-1)*T;               % Time vector

%Behandling af data
f = Fs/2*linspace(0,1,L/2+1);

theRec = d1;
Y = fft(theRec)/L;
mag1 = 2*abs(Y(1:L/2+1));
magdB1 = 20*log10(abs(Y(1:L/2+1)));

theRec = d2;
Y = fft(theRec)/L;
mag2 = 2*abs(Y(1:L/2+1));
magdB2 = 20*log10(abs(Y(1:L/2+1)));

theRec = d3;
Y = fft(theRec)/L;
mag3 = 2*abs(Y(1:L/2+1));
magdB3 = 20*log10(abs(Y(1:L/2+1)));


%Plot af 8 bit
figure(1)
subplot(2,2,1),plot(f,mag1,f,mag2,f,mag3);grid;
title('Single-Sided Amplitude Spectrum')
xlabel('Frequency (Hz)')
ylabel('Magnitude')

subplot(2,2,2),plot(f,magdB1,f,magdB2,f,magdB3);grid;
title('Single-Sided Amplitude Spectrum')
xlabel('Frequency (Hz)')
ylabel('Magnitude in dB')


%Udregning af SNR
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

peak = mag2([PF1 PF2 PF3 PF4 PF5 PF6 PF7 PF8]);
SNR3  = 20*log10((3*peak(1))/(peak(2)+peak(3)+peak(4)));
SNR4  = 20*log10((3*peak(5))/(peak(6)+peak(7)+peak(8)));

peak = mag3([PF1 PF2 PF3 PF4 PF5 PF6 PF7 PF8]);
SNR5  = 20*log10((3*peak(1))/(peak(2)+peak(3)+peak(4)));
SNR6  = 20*log10((3*peak(5))/(peak(6)+peak(7)+peak(8)));

SNR = [SNR1 SNR2; SNR3 SNR4; SNR5 SNR6];

subplot(2,1,2),bar(SNR, 'grouped');grid;
title('SNR')
xlabel('Test')
ylabel('SNR in dB')
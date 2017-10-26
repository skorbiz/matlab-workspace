%Data behandling til Detectoring af en tone tone:
%clear;
% 8 bit - 8000 Hz
% 8 bit - 22050 Hz
% 16 bit - 8000 Hz
% 16 bit - 22050 Hz
% *** FAST BUFFER PÅ 1000 SAMPLES ***
% Ved fast buffer størrelse er det vigtigere at få flere perioder, end en
% høj opløsning



%Extract data
d1 = DetectorTest0x2D8bit0x2D8000';
d2 = DetectorTest0x2D8bit0x2D22050';
d3 = DetectorTest0x2D16bit0x2D8000';
d4 = DetectorTest0x2D16bit0x2D22050';


%Data omkring samplinger
Fs1 = 8000;                    % Sampling frequency
Fs2 = 22050;                   % Sampling frequency
T1 = 1/Fs1;                    % Sample time
T2 = 1/Fs2;                    % Sample time
L = 1000;                      % Length of signal
t1 = (0:L-1)*T1;               % Time vector
t2 = (0:L-1)*T2;               % Time vector

%Behandling af data
f1 = Fs1/2*linspace(0,1,L/2+1);
f2 = Fs2/2*linspace(0,1,L/2+1);

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

theRec = d4;
Y = fft(theRec)/L;
mag4 = 2*abs(Y(1:L/2+1));
magdB4 = 20*log10(abs(Y(1:L/2+1)));


%Plot af 8 bit
figure(1)
subplot(3,2,1),plot(f1,mag1,f2,mag2);grid;axis([0 4100 0 30]);
title('Single-Sided Amplitude Spectrum')
xlabel('Frequency (Hz)')
ylabel('Magnitude')
subplot(3,2,2),plot(f1,magdB1,f2,magdB2);grid;axis([0 4100 -50 50]);
title('Single-Sided Amplitude Spectrum')
xlabel('Frequency (Hz)')
ylabel('Magnitude in dB')

%Plot af 16 bit
subplot(3,2,3),plot(f1,mag3,f2,mag4);grid;axis([0 4100 0 10000]);
title('Single-Sided Amplitude Spectrum')
xlabel('Frequency (Hz)')
ylabel('Magnitude')
subplot(3,2,4),plot(f1,magdB3,f2,magdB4);grid;axis([0 4100 -50 100]);
title('Single-Sided Amplitude Spectrum')
xlabel('Frequency (Hz)')
ylabel('Magnitude in dB')

%Udregning af SNR
PF1_1 = uint64(697*L/Fs1+1);
PF2_1 = uint64(770*L/Fs1+1);
PF3_1 = uint64(852*L/Fs1+1);
PF4_1 = uint64(941*L/Fs1+1);
PF5_1 = uint64(1209*L/Fs1+1);
PF6_1 = uint64(1336*L/Fs1+1);
PF7_1 = uint64(1477*L/Fs1+1);
PF8_1 = uint64(1633*L/Fs1+1);

PF1_2 = uint64(697*L/Fs2+1);
PF2_2 = uint64(770*L/Fs2+1);
PF3_2 = uint64(852*L/Fs2+1);
PF4_2 = uint64(941*L/Fs2+1);
PF5_2 = uint64(1209*L/Fs2+1);
PF6_2 = uint64(1336*L/Fs2+1);
PF7_2 = uint64(1477*L/Fs2+1);
PF8_2 = uint64(1633*L/Fs2+1);


peak = mag1([PF1_1 PF2_1 PF3_1 PF4_1 PF5_1 PF6_1 PF7_1 PF8_1]);
SNR1  = 20*log10((3*peak(1))/(peak(2)+peak(3)+peak(4)));
SNR2  = 20*log10((3*peak(5))/(peak(6)+peak(7)+peak(8)));

peak = mag2([PF1_2 PF2_2 PF3_2 PF4_2 PF5_2 PF6_2 PF7_2 PF8_2]);
SNR3  = 20*log10((3*peak(1))/(peak(2)+peak(3)+peak(4)));
SNR4  = 20*log10((3*peak(5))/(peak(6)+peak(7)+peak(8)));

peak = mag3([PF1_1 PF2_1 PF3_1 PF4_1 PF5_1 PF6_1 PF7_1 PF8_1]);
SNR5  = 20*log10((3*peak(1))/(peak(2)+peak(3)+peak(4)));
SNR6  = 20*log10((3*peak(5))/(peak(6)+peak(7)+peak(8)));

peak = mag4([PF1_2 PF2_2 PF3_2 PF4_2 PF5_2 PF6_2 PF7_2 PF8_2]);
SNR7  = 20*log10((3*peak(1))/(peak(2)+peak(3)+peak(4)));
SNR8  = 20*log10((3*peak(5))/(peak(6)+peak(7)+peak(8)));

SNR = [SNR1 SNR2; SNR3 SNR4; SNR5 SNR6; SNR7 SNR8];

subplot(3,1,3),bar(SNR, 'grouped');grid;
title('SNR')
ylabel('SNR in dB')
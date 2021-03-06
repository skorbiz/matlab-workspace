%Data behandling til generate tone:

%clear;

%%{
%Data omkring samplinger
Fs = 16000;                  % Sampling frequency
T = 1/Fs;                    % Sample time
L = 2000;                    % Length of signal
t = (0:L-1)*T;               % Time vector
nBits = 16;                  % Number of bits used in sample  
nChanels = 1;                % Channels
%%}

%Extract data
d1 = data(:,1);             %Tone 0, 8 bit, 8kHz
d2 = data(:,2);
d3 = data(:,3);             %Tone 0, 8 bit, 16kHz
d4 = data(:,4);
d5 = data(:,5);             %Tone 0, 16 bit, 8kHz
d6 = data(:,6);
d7 = data(:,7);             %Tone 0, 16 bit, 16kHz
d8 = data(:,8);
d9 = data(:,9);             %Tone 0, 16 bit, 4kHz
d10 = data(:,10);
d11 = data(:,11);           %Tone 15, 8 bit, 4kHz

%Behandling af tone 0
f = Fs/2*linspace(0,1,L/2+1);

theRec = d1;
Y = fft(theRec)/L;
mag1 = 2*abs(Y(1:L/2+1));
magdB1 = 20*log10(abs(Y(1:L/2+1)));

theRec = d3;
Y = fft(theRec)/L;
mag3 = 2*abs(Y(1:L/2+1));
magdB3 = 20*log10(abs(Y(1:L/2+1)));

theRec = d5;
Y = fft(theRec)/L;
mag5 = 2*abs(Y(1:L/2+1));
magdB5 = 20*log10(abs(Y(1:L/2+1)));

theRec = d7;
Y = fft(theRec)/L;
mag7 = 2*abs(Y(1:L/2+1));
magdB7 = 20*log10(abs(Y(1:L/2+1)));

theRec = d9;
Y = fft(theRec)/L;
mag9 = 2*abs(Y(1:L/2+1));
magdB9 = 20*log10(abs(Y(1:L/2+1)));

%Plot af tone 0
figure(3)
subplot(3,1,1),plot(f,mag1,f,mag3,f,mag5,f,mag7,f,mag9);grid;
title('Single-Sided Amplitude Spectrum')
xlabel('Frequency (Hz)')
ylabel('Magnitude')
subplot(3,1,2),plot(f,magdB1,f,magdB3,f,magdB5,f,magdB7,f,magdB9);grid;
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

peak = mag3([PF1 PF2 PF3 PF4 PF5 PF6 PF7 PF8]);
SNR3  = 20*log10((3*peak(1))/(peak(2)+peak(3)+peak(4)));
SNR4  = 20*log10((3*peak(5))/(peak(6)+peak(7)+peak(8)));

peak = mag5([PF1 PF2 PF3 PF4 PF5 PF6 PF7 PF8]);
SNR5  = 20*log10((3*peak(1))/(peak(2)+peak(3)+peak(4)));
SNR6  = 20*log10((3*peak(5))/(peak(6)+peak(7)+peak(8)));

peak = mag7([PF1 PF2 PF3 PF4 PF5 PF6 PF7 PF8]);
SNR7  = 20*log10((3*peak(1))/(peak(2)+peak(3)+peak(4)));
SNR8  = 20*log10((3*peak(5))/(peak(6)+peak(7)+peak(8)));

peak = mag9([PF1 PF2 PF3 PF4 PF5 PF6 PF7 PF8]);
SNR9  = 20*log10((3*peak(1))/(peak(2)+peak(3)+peak(4)));
SNR10  = 20*log10((3*peak(5))/(peak(6)+peak(7)+peak(8)));

SNR = [SNR1 SNR2; SNR3 SNR4; SNR5 SNR6; SNR7 SNR8; SNR9 SNR10];

subplot(3,1,3),bar(SNR, 'grouped');grid;
title('SNR')
xlabel('Test nr.')
ylabel('SNR')

%Behandling af tone 15
theRec = d2;
Y = fft(theRec)/L;
mag2 = 2*abs(Y(1:L/2+1));
magdB2 = 20*log10(abs(Y(1:L/2+1)));

theRec = d4;
Y = fft(theRec)/L;
mag4 = 2*abs(Y(1:L/2+1));
magdB4 = 20*log10(abs(Y(1:L/2+1)));

theRec = d6;
Y = fft(theRec)/L;
mag6 = 2*abs(Y(1:L/2+1));
magdB6 = 20*log10(abs(Y(1:L/2+1)));

theRec = d8;
Y = fft(theRec)/L;
mag8 = 2*abs(Y(1:L/2+1));
magdB8 = 20*log10(abs(Y(1:L/2+1)));

theRec = d10;
Y = fft(theRec)/L;
mag10 = 2*abs(Y(1:L/2+1));
magdB10 = 20*log10(abs(Y(1:L/2+1)));

theRec = d11;
Y = fft(theRec)/L;
mag11 = 2*abs(Y(1:L/2+1));
magdB11 = 20*log10(abs(Y(1:L/2+1)));


%Plot af tone 15
figure(2)
subplot(3,1,1),plot(f,mag2,f,mag4,f,mag6,f,mag8,f,mag10,f,mag11);grid;
title('Single-Sided Amplitude Spectrum')
xlabel('Frequency (Hz)')
ylabel('Magnitude')
subplot(3,1,2),plot(f,magdB2,f,magdB4,f,magdB6,f,magdB8,f,magdB10,f,magdB11);grid;
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

peak = mag2([PF1 PF2 PF3 PF4 PF5 PF6 PF7 PF8]);
SNR1  = 20*log10((3*peak(4))/(peak(2)+peak(3)+peak(1)));
SNR2  = 20*log10((3*peak(8))/(peak(6)+peak(7)+peak(5)));

peak = mag4([PF1 PF2 PF3 PF4 PF5 PF6 PF7 PF8]);
SNR3  = 20*log10((3*peak(4))/(peak(2)+peak(3)+peak(1)));
SNR4  = 20*log10((3*peak(8))/(peak(6)+peak(7)+peak(5)));

peak = mag6([PF1 PF2 PF3 PF4 PF5 PF6 PF7 PF8]);
SNR5  = 20*log10((3*peak(4))/(peak(2)+peak(3)+peak(1)));
SNR6  = 20*log10((3*peak(8))/(peak(6)+peak(7)+peak(5)));

peak = mag8([PF1 PF2 PF3 PF4 PF5 PF6 PF7 PF8]);
SNR7  = 20*log10((3*peak(4))/(peak(2)+peak(3)+peak(1)));
SNR8  = 20*log10((3*peak(8))/(peak(6)+peak(7)+peak(5)));

peak = mag10([PF1 PF2 PF3 PF4 PF5 PF6 PF7 PF8]);
SNR9  = 20*log10((3*peak(4))/(peak(2)+peak(3)+peak(1)));
SNR10  = 20*log10((3*peak(8))/(peak(6)+peak(7)+peak(5)));

peak = mag11([PF1 PF2 PF3 PF4 PF5 PF6 PF7 PF8]);
SNR11  = 20*log10((3*peak(4))/(peak(2)+peak(3)+peak(1)));
SNR12  = 20*log10((3*peak(8))/(peak(6)+peak(7)+peak(5)));


SNR = [SNR1 SNR2; SNR3 SNR4; SNR5 SNR6; SNR7 SNR8; SNR9 SNR10; SNR11 SNR12];

subplot(3,1,3),bar(SNR, 'grouped');grid;
title('SNR')
xlabel('Test nr.')
ylabel('SNR')
clc,clear,close All

%Defining a signal
n               = 4096
time            = linspace(-2,2,n);
signal_sine1    = sin(time(1:n/2)*2*pi*2);
signal_sine2    = sin(time(n/2+1:n)*2*pi*20);
signal_sine     = [signal_sine1 signal_sine2];
white_noise     = (rand([1,n])-0.5)*0.4;
signal          = signal_sine+white_noise;

subplot(4,1,1)
plot(time,signal)
title('Signal')
xlabel('time')
ylabel('amplitude')

subplot(4,1,2)
%scale = 2.^(1:8)
scale = 1:500
c = cwt(signal,scale,'db1','absglb');
title('CWT')

subplot(4,1,3)
plot(1:n, c(250,:))
title('Coeficient line for a = 256')
xlabel('time or scale (b)')
ylabel('amplitude')

subplot(4,1,4)
plot(1:n, c(1,:))
title('Coeficient line for a = 2')
xlabel('time or scale (b)')
ylabel('amplitude')
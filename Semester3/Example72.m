clear;
%Example 7.2
%Magnitude and phase plot
figure (1)
[hz w]= freqz([0.187 0.2 0.187],[1],512);
hz = hz.*hamming(512);
phi = 180*unwrap(angle(hz))/pi;
subplot(3,1,1),plot(w,20*log(abs(hz)));grid;
xlabel('Frequency (rad)');
ylabel('Magnitude (dB)');
subplot(3,1,2),plot(w,phi);grid;
xlabel('Frequency (rad)');
ylabel('Phase (degree)');

%w = [0:0.1:2*pi];
mag = abs(0.2+0.3472*cos(w));
subplot(3,1,3),plot(w,mag.*hamming(512));grid;
xlabel('Frequency (rad)');
ylabel('Manitude');

n=[0:1:16]
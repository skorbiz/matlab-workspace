clear;
%Example 72 extendet
Q = 2*pi*800*(1/8000) %Q = 2pi * fc * Ts
N = 91; % The order of the filter | uneaven number
M =(N-1)/2
n = [1:1:M];
h = sin(Q*n)./(pi*n);
h0 = Q/pi;

b = [h(M-n+1) h0 h(n)]

figure (1)
%b = b.*hamming(N)' %Hamming window applied
[hz w]= freqz([b],[1],2000);
w = w.*8000./(2*pi);        % The frequencies is returned to Hz

phi = 180*unwrap(angle(hz))/pi;
subplot(3,1,1),plot(w,20*log(abs(hz)));grid;
xlabel('Frequency (rad)');
ylabel('Magnitude (dB)');
subplot(3,1,2),plot(w,phi);grid;
xlabel('Frequency (rad)');
ylabel('Phase (degree)');

%w = [0:0.1:2*pi];
%mag = abs(0.2+0.3472*cos(w));
%subplot(3,1,3),plot(w,mag);grid;
%xlabel('Frequency (rad)');
%ylabel('Manitude');

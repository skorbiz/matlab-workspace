clear;
Ql = 800*2*pi/8000;
Qh = 1800*2*pi/8000;
N = 57;           % The order of the filter | uneaven number
M =(N-1)/2
n = [1:1:M];
h = sin(Qh*n)./(pi*n)-sin(Ql*n)./(pi*n);
h0 = (Qh-Ql)/pi;

b = [h(M-n+1) h0 h(n)]

%b = b.*kaiser(N)'       %Kaiser window applied
[hz w]= freqz([b],[1],2000);
w = w.*8000./(2*pi);    %w is calculated to be en Hz insted of Rad
phi = 180*unwrap(angle(hz))/pi;
subplot(2,1,1),plot(w,20*log10(abs(hz)));grid;
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
subplot(2,1,2),plot(w,phi);grid;
xlabel('Frequency (Hz)');
ylabel('Phase (degree)');
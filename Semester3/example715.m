clear
%Example 7.15
fs = 8000;
f = [0 0.15 0.25 0.4 0.5 1];        %Edge frequencies
m = [0 0 1 1 0 0];                  %Ideal magnitudes

w = [18 10 18];
b = firpm(57,f,m,w)
format long;
freqz(b,1,512,fs)
axis([0 fs/2 -80 10])
clc,clear,close All
fs = 8000;
Ts = 1/fs;

%Defining a signal
n               = 4096;
time            = linspace(-2,2,n);
signal_sine1    = sin(time(1:n/2)*2*pi*2);
signal_sine2    = sin(time(n/2+1:n)*2*pi*20);
signal_sine     = [signal_sine1 signal_sine2];
white_noise     = (rand([1,n])-0.5)*0.4;
signal          = signal_sine+white_noise;

signallength = n;
wavelet = 'haar';
maxlevel = 9;

[C,L] = wavedec(signal, maxlevel, wavelet);

%Wavelets
A = zeros(maxlevel, signallength);
D = zeros(maxlevel, signallength);

for j=1:maxlevel
    A(j,:) = wrcoef('a',C,L,wavelet,j);
    D(j,:) = wrcoef('d',C,L,wavelet,j);
end
figure('Name',['Signal and contributes from individually reconstructed approximations (Aj) and details (Dj) for analysis with Wavelet ', wavelet])
subplot(maxlevel+1,2,1)
plot(signal,'r-'),ylabel('signal'),axis([0 signallength-1 min(signal) max(signal)]),title('approximations')
subplot(maxlevel+1,2,2)
plot(signal,'r-'),ylabel('signal'),axis([0 signallength-1 min(signal) max(signal)]),title('details')
for j = 1:maxlevel
    subplot(maxlevel+1,2,2*j+1)
    plot(A(j,:),'b-'),ylabel(['A ' num2str(j)]),axis([0 signallength-1 min(A(j,:)) max(A(j,:))])
    subplot(maxlevel+1,2,2*j+2)
    plot(D(j,:),'k-'),ylabel(['D ' num2str(j)]),axis([0 signallength-1 min(D(j,:)) max(D(j,:))])
end


%Coficients
cA = appcoef(C,L,wavelet,maxlevel);
cD = zeros(maxlevel,L(maxlevel+1));
for j=1:maxlevel
    levellength(j) = L(maxlevel+2-j);
    cD(j,1:levellength(j)) = detcoef(C, L, j);
end

% Signal Processing on DWT/FWT coefficients  (cAJ,cDj ---> cAJ_modified,cDj_modified)

cA_mod = cA;
cD_mod = cD;

approxlevel = 5;  %signal/level 1/5, 2/4, 3/3, 4/3, 5/5, 6/3
for j = 1:approxlevel
	cD_mod(j,1:levellength(j)) = zeros(1,levellength(j));
end

% perform classical filtering noise-reduction for comparison
NLP = 50;
LPdelay = NLP/2;
fc = 1000;                       % signal/fc  1/200, 2/500, 3/500, 4/1000, 5/1000, 6/500
hLP = fir1(NLP, fc/(fs/2));
s_mod_filter = filter(hLP,1,signal);
s_mod_filter = [s_mod_filter(LPdelay+1:length(signal)) zeros(1,LPdelay)];
      

% copy from modified cA-vector and cD-matrix to modified [C,L] structure
% (for use in reconstruction 'wrcoef' commands)
C_mod = cA_mod;
for j=maxlevel:-1:1
    C_mod = [C_mod cD_mod(j,1:levellength(j))];
end

% calculate modified signal
s_mod = wrcoef('a',C_mod,L,wavelet,0);
if noise_reduction_method == 6
    s_mod = s_mod_stationary;
end



%%
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
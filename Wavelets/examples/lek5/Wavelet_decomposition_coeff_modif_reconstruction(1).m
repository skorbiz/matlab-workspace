% illustration of basic Wavelet DSP: DWT + modification of coefficients + IDWT

clc;clear;close all;

% ---------------------------------------------------------------------------------------------

N = 1024;
fs = 8000;
Ts = 1/fs;
t = (0:N-1)*Ts;
A1 = 1;
f1 = fs/256;
x = A1*sin(2*pi*f1*t);
sigma = 0.2;
s = x + sigma*randn(1,N);
signallength = length(s);

% ---------------------------------------------------------------------------------------------

% perform DWT/FWT (signal ---> coefficients cAJ,cDj)
wavelet = 'db1';
maxlevel = 5;
[C,L] = wavedec(s, maxlevel, wavelet);
% copy from [C,L] structure to cA-vector and cD-matrix
cA = appcoef(C,L,wavelet,maxlevel);
cD = zeros(maxlevel,L(maxlevel+1));
for j=1:maxlevel
    levellength(j) = L(maxlevel+2-j);
    cD(j,1:levellength(j)) = detcoef(C, L, j);
end

% plot DFT/FWT coefficients at actual length
figure('Name',['Signal and decomposition in DWT/FWT coefficients for analysis with Wavelet ' wavelet , '  (coefficients shown at correct length)'])
subplot(maxlevel+2,2^maxlevel,[1:2^maxlevel])
plot(s,'r-','LineWidth',1.5),ylabel('signal'),axis([0 signallength-1 min(s) max(s)])
for j=1:maxlevel
    subplot(maxlevel+2,2^maxlevel,[j*2^maxlevel+(1:2^(maxlevel-j))])
    plot(cD(j,1:levellength(j)),'k-','LineWidth',1.5),ylabel(['cD ' num2str(j)]),axis([0 levellength(j)-1 min(cD(j,1:levellength(j))) max(cD(j,1:levellength(j)))])
end
subplot(maxlevel+2,2^maxlevel,(maxlevel+1)*2^maxlevel+1)
plot(cA,'b-','LineWidth',1.5),ylabel(['cA ' num2str(maxlevel)]),axis([0 length(cA)-1 min(cA) max(cA)])

% plot DFT/FWT coefficients expanded to signal length
figure('Name',['Signal and decomposition in DWT/FWT coefficients for analysis with Wavelet ' wavelet , '  (coefficients expanded to signal length)'])
subplot(maxlevel+2,1,1)
plot(s,'r-','LineWidth',1.5),ylabel('signal'),axis([0 signallength-1 min(s) max(s)])
for j=1:maxlevel
    subplot(maxlevel+2,1,j+1)
    plot(cD(j,1:levellength(j)),'k-','LineWidth',1.5),ylabel(['cD ' num2str(j)]),axis([0 levellength(j)-1 min(cD(j,1:levellength(j))) max(cD(j,1:levellength(j)))])
end
subplot(maxlevel+2,1,maxlevel+2)
plot(cA,'b-','LineWidth',1.5),ylabel(['cA ' num2str(maxlevel)]),axis([0 length(cA)-1 min(cA) max(cA)])

% perform "hypothetical" IDWT/IFWT (coefficients cAj,cDj ---> Aj,Dj)
% to show contributes to signal from individual details and approximations
A = zeros(maxlevel,signallength);
D = zeros(maxlevel,signallength);
for j=1:maxlevel
    A(j,:) = wrcoef('a',C,L,wavelet,j);
    D(j,:) = wrcoef('d',C,L,wavelet,j);
end

% plot reconstructed approximations and details
figure('Name',['Signal and contributes from individually reconstructed approximations (Aj) and details (Dj) for analysis with Wavelet ', wavelet])
subplot(maxlevel+1,2,1)
plot(s,'r-','LineWidth',1.5),ylabel('signal'),axis([0 signallength-1 min(s) max(s)]),title('approximations')
subplot(maxlevel+1,2,2)
plot(s,'r-','LineWidth',1.5),ylabel('signal'),axis([0 signallength-1 min(s) max(s)]),title('details')
for j = 1:maxlevel
    subplot(maxlevel+1,2,2*j+1)
    plot(A(j,:),'b-','LineWidth',1.5),ylabel(['A ' num2str(j)]),axis([0 signallength-1 min(A(j,:)) max(A(j,:))])
    subplot(maxlevel+1,2,2*j+2)
    plot(D(j,:),'k-','LineWidth',1.5),ylabel(['D ' num2str(j)]),axis([0 signallength-1 min(D(j,:)) max(D(j,:))])
end

% ---------------------------------------------------------------------------------------------


% Signal Processing on DWT/FWT coefficients  (cAJ,cDj ---> cAJ_modified,cDj_modified)
% set cD1 = 0 for second, third and fourth quarter of signal
% set cD2 = 0 for third and fourth quarter of signal
% set cD3 = 0 for fourth quarter of signal

cA_mod = cA;
cD_mod = cD;
cD_mod(1,floor(0.25*levellength(1)):floor(1.00*levellength(1))) = 0;
cD_mod(2,floor(0.50*levellength(2)):floor(1.00*levellength(2))) = 0;
cD_mod(3,floor(0.75*levellength(3)):floor(1.00*levellength(3))) = 0;

% copy from modified cA-vector and cD-matrix to modified [C,L] structure
% (for use in reconstruction 'wrcoef' commands)
C_mod = cA_mod;
for j=maxlevel:-1:1
    C_mod = [C_mod cD_mod(j,1:levellength(j))];
end

% calculate modified signal by IFWT
s_mod = wrcoef('a',C_mod,L,wavelet,0);


% ---------------------------------------------------------------------------------------------

% plot modified DFT/FWT coefficients at actual length
figure('Name','Modified signal and modified DWT/FWT coefficients  (coefficients shown at correct length)')
subplot(maxlevel+2,2^maxlevel,[1:2^maxlevel])
plot(s_mod,'r-','LineWidth',1.5),ylabel('signal{\_}mod'),axis([0 signallength-1 min(s_mod) max(s_mod)])
for j=1:maxlevel
    subplot(maxlevel+2,2^maxlevel,[j*2^maxlevel+(1:2^(maxlevel-j))])
    plot(cD_mod(j,1:levellength(j)),'k-','LineWidth',1.5),ylabel(['cD{\_}mod ' num2str(j)]),axis([0 levellength(j)-1 min(cD_mod(j,1:levellength(j))) max(cD_mod(j,1:levellength(j)))])
end
subplot(maxlevel+2,2^maxlevel,(maxlevel+1)*2^maxlevel+1)
plot(cA_mod,'b-','LineWidth',1.5),ylabel(['cA{\_}mod ' num2str(maxlevel)]),axis([0 length(cA_mod)-1 min(cA_mod) max(cA_mod)])

% plot modified DFT/FWT coefficients expanded to signal length
figure('Name','Modified signal and modified DWT/FWT coefficients  (coefficients expanded to signal length)')
subplot(maxlevel+2,1,1)
plot(s_mod,'r-','LineWidth',1.5),ylabel('signal{\_}mod'),axis([0 signallength-1 min(s_mod) max(s_mod)])
for j=1:maxlevel
    subplot(maxlevel+2,1,j+1)
    plot(cD_mod(j,1:levellength(j)),'k-','LineWidth',1.5),ylabel(['cD{\_}mod ' num2str(j)]),axis([0 levellength(j)-1 min(cD_mod(j,1:levellength(j))) max(cD_mod(j,1:levellength(j)))])
end
subplot(maxlevel+2,1,maxlevel+2)
plot(cA_mod,'b-','LineWidth',1.5),ylabel(['cA{\_}mod ' num2str(maxlevel)]),axis([0 length(cA_mod)-1 min(cA_mod) max(cA_mod)])

% perform "real" IDWT/IFWT (coefficients cAJ_modified,cDj_modified ---> Aj_modified,Dj_modified --->)
% to show contributes to signal from individual modified details and approximations
A_mod = zeros(maxlevel,signallength);
D_mod = zeros(maxlevel,signallength);
for j=1:maxlevel
    A_mod(j,:) = wrcoef('a',C_mod,L,wavelet,j);
    D_mod(j,:) = wrcoef('d',C_mod,L,wavelet,j);
end

% plot modified reconstructed approximations and details
figure('Name',['Modified signal and contributes from modified individually reconstructed approximations (Aj) and details (Dj) for analysis with Wavelet ', wavelet])
subplot(maxlevel+1,2,1)
plot(s_mod,'r-','LineWidth',1.5),ylabel('signal{\_}mod'),axis([0 signallength-1 min(s_mod) max(s_mod)]),title('approximations')
subplot(maxlevel+1,2,2)
plot(s_mod,'r-','LineWidth',1.5),ylabel('signal{\_}mod'),axis([0 signallength-1 min(s_mod) max(s_mod)]),title('details')
for j = 1:maxlevel
    subplot(maxlevel+1,2,2*j+1)
    plot(A_mod(j,:),'b-','LineWidth',1.5),ylabel(['A{\_}mod ' num2str(j)]),axis([0 signallength-1 min(A_mod(j,:)) max(A_mod(j,:))])
    subplot(maxlevel+1,2,2*j+2)
    plot(D_mod(j,:),'k-','LineWidth',1.5),ylabel(['D{\_}mod ' num2str(j)]),axis([0 signallength-1 min(D_mod(j,:)) max(D_mod(j,:))])
end

% plot signal = A0 and modified signal = A0_modified
figure('Name', ['Signal and Wavelet DSP-processed signal     (wavelet = ', wavelet, ')'])
subplot(2,1,1)
plot(s,'r-','LineWidth',1.5),title('signal'),axis([0 signallength-1 min(s) max(s)]), hold on, plot(x,'k-','LineWidth',1.5), hold off
subplot(2,1,2)
plot(s_mod,'b-','LineWidth',1.5),title('modified signal'),axis([0 signallength-1 min(s_mod) max(s_mod)]), hold on, plot(x,'k-','LineWidth',1.5), hold off

% calculate N(0,sigma^2)-equivalent noise densities and SNR before and after processing
sigma_theoretical = sigma
sigma_in = norm(s-x)/sqrt(N)
sigma_out = norm(s_mod-x)/sqrt(N)

SNR_in = (norm(x)/norm(s-x))^2
SNR_out = (norm(x)/norm(s_mod-x))^2

SNR_in_dB = 10*log10(SNR_in)
SNR_out_dB = 10*log10(SNR_out)

% illustration of various Wavelet DSP lossy compression schemes for various signals

clc;clear;close all;
fs = 8000;
Ts = 1/fs;

% ---------------------------------------------------------------------------------------------

signaltype = 5;

switch signaltype
    case {1,'LF-SINE'}
        A1 = 1;
        f1 = 50;
        N = 1024;
        t = (0:N-1)*Ts;
        x = A1*sin(2*pi*f1*t);  
      
    case {2,'SUM-OF-SINES'}
        A1 = 1;
        f1 = 50;
        A2 = 0.2;
        f2 = 250;
        N = 1024;
        t = (0:N-1)*Ts;
        x = A1*sin(2*pi*f1*t)+A2*sin(2*pi*f2*t);
      
    case {3,'CHIRP-SIGNAL'}
        A = 1;
        f1 = 10;
        f2 = 300;
        N = 1024;
        t = (0:N-1)*Ts;
        fase = 2*pi*(f1*t + 0.5*(f2-f1)*t.^2/(N*Ts));
        x = sin(fase);  
       
    case {4,'HF-SPECTRUM'}
        N = 1024;
        data = load('noisbump.mat');
        x_data = data.noisbump;
        hLP = fir1(50,0.1);
        x = filter(hLP,1,x_data);
        x = x_data;
        
    case {5,'MULTILEVEL DIGITAL BASEBAND SIGNAL'}
        N = 1024;
        Nsymbol = 32;
        data = floor(4*rand(1,Nsymbol))-3/2;
        repeat_row = [ones(1,32) zeros(1,1024-32)];
        repeat_matrix(1,:) = repeat_row;
        for p=2:Nsymbol
            repeat_matrix(p,:)=circshift(repeat_row,(p-1)*[0 32]);
        end
        x = data*repeat_matrix;
        
    case {6,'PULSESHAPED MULTILEVEL DIGITAL BASEBAND SIGNAL'}
        N = 1024;
        Nsymbol = 32;
        data = floor(4*rand(1,Nsymbol))-3/2;
        repeat_row = [gausswin(32,2.5)' zeros(1,1024-32)];
        repeat_matrix(1,:) = repeat_row;
        for p=2:Nsymbol
            repeat_matrix(p,:)=circshift(repeat_row,(p-1)*[0 32]);
        end
        x = data*repeat_matrix;
        
    case {7,'VHF-SPECTRUM'}
        N = 1024;
        data = load('noisbump.mat');
        x_data = data.noisbump;
        hLP = fir1(50,0.1);
        x = filter(hLP,1,x_data);
        x = x_data;
        xd = downsample(x,4);
        x = [xd xd xd xd];  
        
    case {8,'UHF-SPECTRUM'}   
        N = 1024;
        data = load('noisbump.mat');
        x_data = data.noisbump;
        hHP = fir1(50,0.05,'high');
        x = filter(hHP,1,x_data);
             
    case {9,'BARCODE-SIGNAL'}   % a la 1D fingeraftryk: rand evt. et antal samples pr streg
        N = 1024;
        x = [];
        for n = 1:floor(N/4)
            temp = rand;
            bit = (temp >= 0.95);
            x = [x bit bit bit bit];
        end
        
end

signallength = length(x);

% ---------------------------------------------------------------------------------------------

noisetype = 0;     % not relevant for compression studies, see m-file concerning noise-reduction
switch noisetype
    case {0,'NONE'}
      noise = zeros(1,signallength);  
    case {1,'AWGN'}
      sigma = 0.5;
      noise = sigma*randn(1,signallength);  
      
    case {2,'UNIFORM'}
    case {3,'LPF-AWGN'}
    case {4,'NON-STATIONARY AWGN'}
end

% ---------------------------------------------------------------------------------------------

s = x + noise;

% ---------------------------------------------------------------------------------------------

% perform DWT/FWT (signal ---> coefficients cAJ,cDj)
wavelet = 'db1';    %signal/wavelet  1/db6, 2/db4, 3/db10, 4/db5, 5/db1, 6/db6, 7/db3, 8/db2, 9/db1
compresslevel = 2;  %signal/level 1/5, 2/4, 3/3, 4/2, 5/5, 6/3, 7/5, 8/2, 9/2
[C,L] = wavedec(s, compresslevel, wavelet);
% copy from [C,L] structure to cA-vector and cD-matrix
cA = appcoef(C,L,wavelet,compresslevel);
cD = zeros(compresslevel,L(compresslevel+1));
for j=1:compresslevel
    levellength(j) = L(compresslevel+2-j);
    cD(j,1:levellength(j)) = detcoef(C, L, j);
end

% plot DWT/FWT coefficients at actual length
figure('Name',['Signal and decomposition in DWT/FWT coefficients for analysis with Wavelet ' wavelet , '  (coefficients shown at correct length)'])
subplot(compresslevel+2,2^compresslevel,[1:2^compresslevel])
plot(s,'r-'),ylabel('signal'),axis([0 signallength-1 min(s) max(s)])
for j=1:compresslevel
    subplot(compresslevel+2,2^compresslevel,[j*2^compresslevel+(1:2^(compresslevel-j))])
    plot(cD(j,1:levellength(j)),'k-'),ylabel(['cD ' num2str(j)]),xlim([0 levellength(j)-1])
end
subplot(compresslevel+2,2^compresslevel,(compresslevel+1)*2^compresslevel+1)
plot(cA,'b-'),ylabel(['cA ' num2str(compresslevel)]),axis([0 length(cA)-1 min(cA) max(cA)])

% plot DWT/FWT coefficients expanded to signal length
figure('Name',['Signal and decomposition in DWT/FWT coefficients for analysis with Wavelet ' wavelet , '  (coefficients expanded to signal length)'])
subplot(compresslevel+2,1,1)
plot(s,'r-'),ylabel('signal'),axis([0 signallength-1 min(s) max(s)])
for j=1:compresslevel
    subplot(compresslevel+2,1,j+1)
    plot(cD(j,1:levellength(j)),'k-'),ylabel(['cD ' num2str(j)]),xlim([0 levellength(j)-1])
end
subplot(compresslevel+2,1,compresslevel+2)
plot(cA,'b-'),ylabel(['cA ' num2str(compresslevel)]),axis([0 length(cA)-1 min(cA) max(cA)])

% ---------------------------------------------------------------------------------------------

% perform "hypothetical" IDWT/IFWT (coefficients cAj,cDj ---> Aj,Dj)
% to show contributes to signal from individual details and approximations
A = zeros(compresslevel,signallength);
D = zeros(compresslevel,signallength);
for j=1:compresslevel
    A(j,:) = wrcoef('a',C,L,wavelet,j);
    D(j,:) = wrcoef('d',C,L,wavelet,j);
end

% plot reconstructed approximations and details
figure('Name',['Signal and contributes from individually reconstructed approximations (Aj) and details (Dj) for analysis with Wavelet ', wavelet])
subplot(compresslevel+1,2,1)
plot(s,'r-'),ylabel('signal'),axis([0 signallength-1 min(s) max(s)]),title('approximations')
subplot(compresslevel+1,2,2)
plot(s,'r-'),ylabel('signal'),axis([0 signallength-1 min(s) max(s)]),title('details')
for j = 1:compresslevel
    subplot(compresslevel+1,2,2*j+1)
    plot(A(j,:),'b-'),ylabel(['A ' num2str(j)]),axis([0 signallength-1 min(A(j,:)) max(A(j,:))])
    subplot(compresslevel+1,2,2*j+2)
    plot(D(j,:),'k-'),ylabel(['D ' num2str(j)]),xlim([0 signallength-1])
end

% ---------------------------------------------------------------------------------------------
% ---------------------------------------------------------------------------------------------
% ---------------------------------------------------------------------------------------------

% Signal Processing (compression) on DWT/FWT coefficients  (cAJ,cDj ---> cAJ_modified,cDj_modified)

cA_mod = cA;
cD_mod = cD;

compression_method = 1;
switch compression_method
    case {0,'NONE'}
                
    case {1,'SELECT APPROXIMATION = REMOVE HF-DETAILS'}
        for j = 1:compresslevel
            cD_mod(j,1:levellength(j)) = zeros(1,levellength(j));
        end
        
    case {2,'GLOBAL HARD THRESHOLDING OF DETAIL COEFFS'}
        switch signaltype
            case 1                                        % supply one global threshold for level 1 to compresslevel
                threshold = 0.4483;
            case 2
                threshold = 0.3867;
            case 3
                threshold = 0.3829;
            case 4
                threshold = 2.695;  % level 3
            case 5
                threshold = 0;
            case 6
                threshold = 0.376;
            case 7
                threshold = 3.241;
            case 9
                threshold = 0.3;    % level 5
        end
        for j = 1:compresslevel
            lower_cD = find(abs(cD(j,1:levellength(j))) <= threshold);
            cD_mod(j,lower_cD) = 0;
        end
        
    case {3,'LEVELWISE GLOBAL HARD THRESHOLDING OF DETAIL COEFFS'}
        switch signaltype
            case 1                                        % supply thresholds from level 1 to compresslevel
                thresholds = [0.005 0.005 0.034 0.046 0.074];
            case 2
                thresholds = [0.005 0.081 0.044 0.044];
            case 3
                thresholds = [0.012 0.049 0.061];
            case 4
                thresholds = [0.0037 0.076 2.346];
            case 5
                thresholds = [0 0 0 0 0 ];
            case 6
                thresholds = [2 2 2 0 0 0 0 0 0 ];
            case 7
                thresholds = [7 7 4 2 1];
        end
        for j = 1:compresslevel
            lower = find(abs(cD(j,1:levellength(j))) <= thresholds(j));
            cD_mod(j,lower) = 0;
        end
end

% copy from modified cA-vector and cD-matrix to modified [C,L] structure
% (for use in reconstruction 'wrcoef' commands)
C_mod = cA_mod;
for j=compresslevel:-1:1
    C_mod = [C_mod cD_mod(j,1:levellength(j))];
end

% calculate reconstructed compressed signal
s_mod = wrcoef('a',C_mod,L,wavelet,0);


% perform classical DFT-based lossy compression for comparison
S_FFT = fft(s,N)/sqrt(N);
S_mod_FFT = S_FFT;
FFT_threshold = 0.1;   % signal/threshold 1/0.05, 4/1, 5/0.5, 9/0.1
lower = find(abs(S_FFT) < FFT_threshold);
S_mod_FFT(lower) = 0;
s_mod_FFT = ifft(S_mod_FFT,N)*sqrt(N);
figure('Name','Classical DFT/FFT and DCT-based compression')
subplot(2,2,1)
plot(abs(S_FFT),'r-'),title('Signal DFT/FFT amplitude spectrum'),xlim([0 N-1])
subplot(2,2,3)
plot(abs(S_mod_FFT),'b-'),title('Compressed signal DFT/FFT amplitude spectrum'),xlim([0 N-1])


% perform classical DCT-based lossy compression for comparison
S_DCT = dct(s,N);
S_mod_DCT = S_DCT;
DCT_threshold = 0.1;   % signal/threshold 1/0.05, 4/1, 5/0.5, 9/0.1
lower = find(abs(S_DCT) < DCT_threshold);
S_mod_DCT(lower) = 0;
s_mod_DCT = idct(S_mod_DCT,N);
subplot(2,2,2)
plot(abs(S_DCT),'r-'),title('Signal DCT amplitude spectrum'),xlim([0 N-1])
subplot(2,2,4)
plot(abs(S_mod_DCT),'b-'),title('Compressed signal DCT amplitude spectrum'),xlim([0 N-1])


% ---------------------------------------------------------------------------------------------
% ---------------------------------------------------------------------------------------------
% ---------------------------------------------------------------------------------------------

% plot modified DWT/FWT coefficients at actual length
figure('Name','Modified signal and modified DWT/FWT coefficients  (coefficients shown at correct length)')
subplot(compresslevel+2,2^compresslevel,[1:2^compresslevel])
plot(s_mod,'r-'),ylabel('signal{\_}mod'),axis([0 signallength-1 min(s_mod) max(s_mod)])
for j=1:compresslevel
    subplot(compresslevel+2,2^compresslevel,[j*2^compresslevel+(1:2^(compresslevel-j))])
    plot(cD_mod(j,1:levellength(j)),'k-'),ylabel(['cD{\_}mod ' num2str(j)]),xlim([0 levellength(j)-1])
end
subplot(compresslevel+2,2^compresslevel,(compresslevel+1)*2^compresslevel+1)
plot(cA_mod,'b-'),ylabel(['cA{\_}mod ' num2str(compresslevel)]),axis([0 length(cA_mod)-1 min(cA_mod) max(cA_mod)])

% plot modified DWT/FWT coefficients expanded to signal length
figure('Name','Modified signal and modified DWT/FWT coefficients  (coefficients expanded to signal length)')
subplot(compresslevel+2,1,1)
plot(s_mod,'r-'),ylabel('signal{\_}mod'),axis([0 signallength-1 min(s_mod) max(s_mod)])
for j=1:compresslevel
    subplot(compresslevel+2,1,j+1)
    plot(cD_mod(j,1:levellength(j)),'k-'),ylabel(['cD{\_}mod ' num2str(j)]),xlim([0 levellength(j)-1])
end
subplot(compresslevel+2,1,compresslevel+2)
plot(cA_mod,'b-'),ylabel(['cA{\_}mod ' num2str(compresslevel)]),axis([0 length(cA_mod)-1 min(cA_mod) max(cA_mod)])

% ---------------------------------------------------------------------------------------------

% perform "real" IDWT/IFWT (coefficients cAJ_modified,cDj_modified ---> Aj_modified,Dj_modified --->)
% to show contributes to signal from individual modified details and approximations
A_mod = zeros(compresslevel,signallength);
D_mod = zeros(compresslevel,signallength);
for j=1:compresslevel
    A_mod(j,:) = wrcoef('a',C_mod,L,wavelet,j);
    D_mod(j,:) = wrcoef('d',C_mod,L,wavelet,j);
end

% plot modified reconstructed approximations and details
figure('Name',['Modified signal and contributes from modified individually reconstructed approximations (Aj) and details (Dj) for analysis with Wavelet ', wavelet])
subplot(compresslevel+1,2,1)
plot(s_mod,'r-'),ylabel('signal{\_}mod'),axis([0 signallength-1 min(s_mod) max(s_mod)]),title('approximations')
subplot(compresslevel+1,2,2)
plot(s_mod,'r-'),ylabel('signal{\_}mod'),axis([0 signallength-1 min(s_mod) max(s_mod)]),title('details')
for j = 1:compresslevel
    subplot(compresslevel+1,2,2*j+1)
    plot(A_mod(j,:),'b-'),ylabel(['A{\_}mod ' num2str(j)]),axis([0 signallength-1 min(A_mod(j,:)) max(A_mod(j,:))])
    subplot(compresslevel+1,2,2*j+2)
    plot(D_mod(j,:),'k-'),ylabel(['D{\_}mod ' num2str(j)]),xlim([0 signallength-1])
end

% ---------------------------------------------------------------------------------------------

% calculate retained energy and compression-degree for current decomposition-compression-reconstruction scheme

energy_coeff_s = (norm(cA))^2;
for j=1:compresslevel
    energy_coeff_s = energy_coeff_s + (norm(cD(j,1:levellength(j))))^2;
end
energy_coeff_s_mod = (norm(cA_mod))^2;
for j=1:compresslevel
    energy_coeff_s_mod = energy_coeff_s_mod + (norm(cD_mod(j,1:levellength(j))))^2;
end
retained_energy_coeff = 100*energy_coeff_s_mod/energy_coeff_s;

energy_coeff_s_FFT = norm(S_FFT)^2;
energy_coeff_s_mod_FFT = norm(S_mod_FFT)^2;
retained_energy_coeff_FFT = 100*energy_coeff_s_mod_FFT/energy_coeff_s_FFT;

energy_coeff_s_DCT = norm(S_DCT)^2;
energy_coeff_s_mod_DCT = norm(S_mod_DCT)^2;
retained_energy_coeff_DCT = 100*energy_coeff_s_mod_DCT/energy_coeff_s_DCT;

relative_L2_norm_error = 100*(norm(s_mod - s)/norm(s))^2;
relative_L2_norm_error_FFT = 100*(norm(s_mod_FFT - s)/norm(s))^2;
relative_L2_norm_error_DCT = 100*(norm(s_mod_DCT - s)/norm(s))^2;

coeff_s = sum(levellength)+levellength(compresslevel);     
zero_coeff_s_mod = 0;
for j = 1:compresslevel
    zero_coeff_s_mod = zero_coeff_s_mod + length(find(cD_mod(j,1:levellength(j)) == 0));
end
zero_coeff = 100*zero_coeff_s_mod/coeff_s;
compressfactor = coeff_s/(coeff_s - zero_coeff_s_mod);

coeff_FFT = N;
zero_coeff_s_mod_FFT = length(find(S_mod_FFT == 0));
zero_coeff_FFT = 100*zero_coeff_s_mod_FFT/coeff_FFT;
compressfactor_FFT = coeff_FFT/(coeff_FFT - zero_coeff_s_mod_FFT);

coeff_DCT = N;
zero_coeff_s_mod_DCT = length(find(S_mod_DCT == 0));
zero_coeff_DCT = 100*zero_coeff_s_mod_DCT/coeff_DCT;
compressfactor_DCT = coeff_DCT/(coeff_DCT - zero_coeff_s_mod_DCT);

% ---------------------------------------------------------------------------------------------

% plot signal and compressed signals
figure('Name', ['   Original signal and compressed signals     (wavelet = ' wavelet ', compresslevel = ' num2str(compresslevel) ')'])     
subplot(4,1,1)
plot(s,'k-','LineWidth',2),title('original signal'),axis([0 signallength-1 min(s) max(s)])
subplot(4,1,2)
plot(s,'k-'),title(['original and WAVELET DSP compressed signal,   (Retained coeff.energy = ' num2str(retained_energy_coeff,4)  '%,  Relative signal L2-norm error = ' num2str(relative_L2_norm_error,4) '%,  Zero-coeffs = '  num2str(zero_coeff,4) '%,  Compress-factor = ' num2str(compressfactor,4) ')']),axis([0 signallength-1 min(s) max(s)]), hold on, plot(s_mod,'r-','LineWidth',2), hold off,legend('original signal','compressed signal')
subplot(4,1,3)
plot(s,'k-'),title(['original and CLASSICAL DFT compressed signal,   (Retained coeff.energy = ' num2str(retained_energy_coeff_FFT,4)  '%,  Relative signal L2-norm error = ' num2str(relative_L2_norm_error_FFT,4) '%,  Zero-coeffs = '  num2str(zero_coeff_FFT,4) '%,  Compress-factor = ' num2str(compressfactor_FFT,4) ')']),axis([0 signallength-1 min(s) max(s)]), hold on, plot(s_mod_FFT,'b-','LineWidth',2), hold off,legend('original signal','compressed signal')
subplot(4,1,4)
plot(s,'k-'),title(['original and CLASSICAL DCT compressed signal,   (Retained coeff.energy = ' num2str(retained_energy_coeff_DCT,4)  '%,  Relative signal L2-norm error = ' num2str(relative_L2_norm_error_DCT,4) '%,  Zero-coeffs = '  num2str(zero_coeff_DCT,4) '%,  Compress-factor = ' num2str(compressfactor_DCT,4) ')']),axis([0 signallength-1 min(s) max(s)]), hold on, plot(s_mod_DCT,'m-','LineWidth',2), hold off,legend('original signal','compressed signal')


% ---------------------------------------------------------------------------------------------

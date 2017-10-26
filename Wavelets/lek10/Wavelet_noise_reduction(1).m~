% illustration of various Wavelet DSP noise reduction schemes for various signals and noise-types

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
        data = load('noisbump.mat');
        x_data = data.noisbump;
        hLP = fir1(50,0.1);
        x = filter(hLP,1,x_data);
        
    case {5,'MULTILEVEL DIGITAL BASEBAND SIGNAL'}
        Nsymbol = 32;
        data = floor(4*rand(1,Nsymbol))-3/2;
        repeat_row = [ones(1,32) zeros(1,1024-32)];
        repeat_matrix(1,:) = repeat_row;
        for p=2:Nsymbol
            repeat_matrix(p,:)=circshift(repeat_row,(p-1)*[0 32]);
        end
        x = data*repeat_matrix;
        
    case {6,'PULSESHAPED MULTILEVEL DIGITAL BASEBAND SIGNAL'}
        Nsymbol = 32;
        data = floor(4*rand(1,Nsymbol))-3/2;
        repeat_row = [gausswin(32,2.5)' zeros(1,1024-32)];
        repeat_matrix(1,:) = repeat_row;
        for p=2:Nsymbol
            repeat_matrix(p,:)=circshift(repeat_row,(p-1)*[0 32]);
        end
        x = data*repeat_matrix;
end

signallength = length(x);

% ---------------------------------------------------------------------------------------------

noisetype = 1;
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
wavelet = 'db1';   %signal/wavelet  1/db6, 2/db4, 3/db10, 4/db5, 5/db1, 6/db6
maxlevel = 9;
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
plot(s,'r-'),ylabel('signal'),axis([0 signallength-1 min(s) max(s)])
for j=1:maxlevel
    subplot(maxlevel+2,2^maxlevel,[j*2^maxlevel+(1:2^(maxlevel-j))])
    plot(cD(j,1:levellength(j)),'k-'),ylabel(['cD ' num2str(j)]),axis([0 levellength(j)-1 min(cD(j,1:levellength(j))) max(cD(j,1:levellength(j)))])
end
subplot(maxlevel+2,2^maxlevel,(maxlevel+1)*2^maxlevel+1)
plot(cA,'b-'),ylabel(['cA ' num2str(maxlevel)]),axis([0 length(cA)-1 min(cA) max(cA)])

% plot DFT/FWT coefficients expanded to signal length
figure('Name',['Signal and decomposition in DWT/FWT coefficients for analysis with Wavelet ' wavelet , '  (coefficients expanded to signal length)'])
subplot(maxlevel+2,1,1)
plot(s,'r-'),ylabel('signal'),axis([0 signallength-1 min(s) max(s)])
for j=1:maxlevel
    subplot(maxlevel+2,1,j+1)
    plot(cD(j,1:levellength(j)),'k-'),ylabel(['cD ' num2str(j)]),axis([0 levellength(j)-1 min(cD(j,1:levellength(j))) max(cD(j,1:levellength(j)))])
end
subplot(maxlevel+2,1,maxlevel+2)
plot(cA,'b-'),ylabel(['cA ' num2str(maxlevel)]),axis([0 length(cA)-1 min(cA) max(cA)])

% ---------------------------------------------------------------------------------------------

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
plot(s,'r-'),ylabel('signal'),axis([0 signallength-1 min(s) max(s)]),title('approximations')
subplot(maxlevel+1,2,2)
plot(s,'r-'),ylabel('signal'),axis([0 signallength-1 min(s) max(s)]),title('details')
for j = 1:maxlevel
    subplot(maxlevel+1,2,2*j+1)
    plot(A(j,:),'b-'),ylabel(['A ' num2str(j)]),axis([0 signallength-1 min(A(j,:)) max(A(j,:))])
    subplot(maxlevel+1,2,2*j+2)
    plot(D(j,:),'k-'),ylabel(['D ' num2str(j)]),axis([0 signallength-1 min(D(j,:)) max(D(j,:))])
end

% ---------------------------------------------------------------------------------------------
% ---------------------------------------------------------------------------------------------
% ---------------------------------------------------------------------------------------------

% Signal Processing on DWT/FWT coefficients  (cAJ,cDj ---> cAJ_modified,cDj_modified)

cA_mod = cA;
cD_mod = cD;

noise_reduction_method = 1;
switch noise_reduction_method
    case {0,'NONE'}
        
    case {1,'SELECT APPROXIMATION = REMOVE HF-DETAILS'}
        approxlevel = 5;  %signal/level 1/5, 2/4, 3/3, 4/3, 5/5, 6/3
        for j = 1:approxlevel
            cD_mod(j,1:levellength(j)) = zeros(1,levellength(j));
        end
        
    case {2,'LEVELWISE GLOBAL HARD THRESHOLDING'}
        switch signaltype
            case 1                                        % supply thresholds from level 1 to maxlevel
                thresholds = [2 1 1 1 1 0 0 0 0];
            case 2
                thresholds = [1 1 1 1 0 0 0 0 0 ];
            case 3
                thresholds = [2 2 2 0 0 0 0 0 0];
            case 4
                thresholds = [3.7 3.7 3.7 3.7 0 0 0 0 0 ];
            case 5
                thresholds = [4 2 2 2 2 0 0 0 0];
            case 6
                thresholds = [2 2 2 0 0 0 0 0 0 ];
        end
        for j = 1:maxlevel
            lower = find(abs(cD(j,1:levellength(j))) <= thresholds(j));
            cD_mod(j,lower) = 0;
        end
        
    case {3,'LEVELWISE GLOBAL SOFT THRESHOLDING'}
        switch signaltype
            case 1                                        % supply thresholds from level 1 to maxlevel
                thresholds = [2 1 1 1 1 0 0 0 0];
            case 2
                thresholds = [1 1 1 1 0 0 0 0 0 ];
            case 3
                thresholds = [2 2 2 0 0 0 0 0 0];
            case 4
                thresholds = [5 5 5 5 0 0 0 0 0];
            case 5
                thresholds = [4 2 2 2 2 0 0 0 0 ];
            case 6
                thresholds = [2 2 2 0 0 0 0 0 0 ];
        end
        for j = 1:maxlevel
            lower = find(abs(cD(j,1:levellength(j))) <= thresholds(j));
            cD_mod(j,lower) = 0;
        end
        for j = 1:maxlevel
            higher_pos = find(cD(j,1:levellength(j)) > thresholds(j));
            cD_mod(j,higher_pos) = cD(j,higher_pos) - thresholds(j);
            higher_neg = find(cD(j,1:levellength(j)) < -thresholds(j));
            cD_mod(j,higher_neg) = cD(j,higher_neg) + thresholds(j);            
        end
              
    case {4,'LOCAL HARD THRESHOLDING'}
    case {5,'LOCAL SOFT THRESHOLDING'}  
        
    case {6,'STATIONARY DISCRETE WAVELET TRANSFORM'}  % only initial studies performed !
        stationary_maxlevel = 3;      % 1/5, 2/3, 3/3, 4/3, 5/2, 6/3
        [swa,swd] = swt(s,stationary_maxlevel,wavelet);
        [thr,sorh] = ddencmp('den','wv',s);
        dswd = wthresh(swd,sorh,thr);
        s_mod_stationary = iswt(swa,dswd,wavelet);       
end

% perform classical filtering noise-reduction for comparison
NLP = 50;
LPdelay = NLP/2;
fc = 1000;                       % signal/fc  1/200, 2/500, 3/500, 4/1000, 5/1000, 6/500
hLP = fir1(NLP, fc/(fs/2));
s_mod_filter = filter(hLP,1,s);
s_mod_filter = [s_mod_filter(LPdelay+1:length(s)) zeros(1,LPdelay)];
      

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

% ---------------------------------------------------------------------------------------------
% ---------------------------------------------------------------------------------------------
% ---------------------------------------------------------------------------------------------

% plot modified DFT/FWT coefficients at actual length
figure('Name','Modified signal and modified DWT/FWT coefficients  (coefficients shown at correct length)')
subplot(maxlevel+2,2^maxlevel,[1:2^maxlevel])
plot(s_mod,'r-'),ylabel('signal{\_}mod'),axis([0 signallength-1 min(s_mod) max(s_mod)])
for j=1:maxlevel
    subplot(maxlevel+2,2^maxlevel,[j*2^maxlevel+(1:2^(maxlevel-j))])
    plot(cD_mod(j,1:levellength(j)),'k-'),ylabel(['cD{\_}mod ' num2str(j)]),axis([0 levellength(j)-1 min(cD(j,1:levellength(j))) max(cD(j,1:levellength(j)))])
end
subplot(maxlevel+2,2^maxlevel,(maxlevel+1)*2^maxlevel+1)
plot(cA_mod,'b-'),ylabel(['cA{\_}mod ' num2str(maxlevel)]),axis([0 length(cA_mod)-1 min(cA_mod) max(cA_mod)])

% plot modified DFT/FWT coefficients expanded to signal length
figure('Name','Modified signal and modified DWT/FWT coefficients  (coefficients expanded to signal length)')
subplot(maxlevel+2,1,1)
plot(s_mod,'r-'),ylabel('signal{\_}mod'),axis([0 signallength-1 min(s_mod) max(s_mod)])
for j=1:maxlevel
    subplot(maxlevel+2,1,j+1)
    plot(cD_mod(j,1:levellength(j)),'k-'),ylabel(['cD{\_}mod ' num2str(j)]),axis([0 levellength(j)-1 min(cD(j,1:levellength(j))) max(cD(j,1:levellength(j)))])
end
subplot(maxlevel+2,1,maxlevel+2)
plot(cA_mod,'b-'),ylabel(['cA{\_}mod ' num2str(maxlevel)]),axis([0 length(cA_mod)-1 min(cA_mod) max(cA_mod)])

% ---------------------------------------------------------------------------------------------

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
plot(s_mod,'r-'),ylabel('signal{\_}mod'),axis([0 signallength-1 min(s_mod) max(s_mod)]),title('approximations')
subplot(maxlevel+1,2,2)
plot(s_mod,'r-'),ylabel('signal{\_}mod'),axis([0 signallength-1 min(s_mod) max(s_mod)]),title('details')
for j = 1:maxlevel
    subplot(maxlevel+1,2,2*j+1)
    plot(A_mod(j,:),'b-'),ylabel(['A{\_}mod ' num2str(j)]),axis([0 signallength-1 min(A_mod(j,:)) max(A_mod(j,:))])
    subplot(maxlevel+1,2,2*j+2)
    plot(D_mod(j,:),'k-'),ylabel(['D{\_}mod ' num2str(j)]),axis([0 signallength-1 min(D(j,:)) max(D(j,:))])
end

% ---------------------------------------------------------------------------------------------

% calculate SNR before and after processing
SNR_in = (norm(x)/norm(s-x))^2;
SNR_out = (norm(x)/norm(s_mod-x))^2;
SNR_out_filter = (norm(x(1:length(s)-LPdelay))/norm(s_mod_filter(1:length(s)-LPdelay)-x(1:length(s)-LPdelay)))^2;
SNR_gain = SNR_out/SNR_in;
SNR_gain_filter = SNR_out_filter/SNR_in;
SNR_in_dB = 10*log10(SNR_in);
SNR_out_dB = 10*log10(SNR_out);
SNR_out_filter_dB = 10*log10(SNR_out_filter);
SNR_gain_dB = SNR_out_dB - SNR_in_dB;
SNR_gain_filter_dB = SNR_out_filter_dB - SNR_in_dB;



% ---------------------------------------------------------------------------------------------

% plot signal = A0 and modified signal = A0_modified
figure('Name', ['Signal and Wavelet DSP-processed signal     (wavelet = ', wavelet, ')'])     
subplot(3,1,1)
plot(x,'k-','LineWidth',4),title(['original signal + noise,  SNR = ' num2str(SNR_in,3)  ' = ' num2str(SNR_in_dB,3) ' dB                                                                                              ']),axis([0 signallength-1 min(s) max(s)]), hold on, plot(s,'r-','LineWidth',1.5), hold off,legend('original signal','original signal + noise')
subplot(3,1,2)
plot(x,'k-'),title(['original and WAVELET DSP estimated signal,  SNR = ' num2str(SNR_out,3)  ' = ' num2str(SNR_out_dB,3) ' dB' '                      (SNR-gain = ' num2str(SNR_gain,3) ' = ' num2str(SNR_gain_dB,3) ' dB)']),axis([0 signallength-1 min(s_mod) max(s_mod)]), hold on, plot(s_mod,'b-','LineWidth',2), hold off,legend('original signal','estimated signal')
subplot(3,1,3)
plot(x,'k-'),title(['original and CLASSICAL DSP estimated signal,  SNR = ' num2str(SNR_out_filter,3)  ' = ' num2str(SNR_out_filter_dB,3) ' dB' '                      (SNR-gain = ' num2str(SNR_gain_filter,3) ' = ' num2str(SNR_gain_filter_dB,3) ' dB)']),axis([0 signallength-1 min(s_mod) max(s_mod)]), hold on, plot(s_mod_filter,'b-','LineWidth',2), hold off,legend('original signal','estimated signal')

% ---------------------------------------------------------------------------------------------

% Illustration of Wavelet analysis capabilities for 
% - detection of various discontinuities in a signal or its derivatives
% - suppresion of polynomial parts of a signal 

% These capabilities are closely linked to two properties of the choseen wavelet for analysis
% - the mathematical "regularity" of the wavelet function
% - the number of vanishing moments of the wavelet function

clc;clear;close all;
fs = 8000;
Ts = 1/fs;
N = 1024;
n = 0:N-1;

% ---------------------------------------------------------------------------------------------

signaltype = 10;
switch signaltype
    
    case {1,'step signal = zero-order derivative discontinuity'}
        x = [zeros(1,500) 0.5 ones(1,N-501)];
        signaldegree = '0';
        signalregularity = '0';
        wavelet = 'db1';
        maxlevel = 6;
        detectlevel = 1;
             
    case {2,'piecewise linear signal with zero-order derivative discontinuity'}    
        x = [400*(500:-1:1)/500 0 10+((1:(N-501))/5)];
        signaldegree = '1';
        signalregularity = '0';
        wavelet = 'db2';
        maxlevel = 6;
        detectlevel = 1;
         
    case {3,'piecewise linear signal with first-order derivative discontinuity'}
        x = [400*(500:-1:1)/500 0 400*(1:(N-501))/500];
        signaldegree = '1';
        signalregularity = 'approx. ---> 1-';
        wavelet = 'db2';
        maxlevel = 6;
        detectlevel = 1;             
        
    case {4,'piecewise parabolic signal with zero-order derivative discontinuity'}
        t = (0:(N/16-1))*Ts;
        x = t.^2;
        x = [x (0.2e-5)+fliplr(x)];
        x = [x x x x x x x x];
        signaldegree = '2';
        signalregularity = '0';
        wavelet = 'db3';
        maxlevel = 6;
        detectlevel = 1;
        
    case {5,'piecewise parabolic signal with first-order derivative discontinuity'}
        t = (0:(N/16-1))*Ts;
        x = t.^2;
        x = [x fliplr(x)];
        x = [x x x x x x x x];
        signaldegree = '2';
        signalregularity = '0.5';
        wavelet = 'db3';
        maxlevel = 6;
        detectlevel = 1;
        
    case {6,'piecewise parabolic signal with second-order derivative discontinuity'}
        t1 = (0:100)*Ts;
        x1 = t1.^2;
        x1 = fliplr(x1);
        t2 = (0:922)*Ts;
        x2 = 0.01*(t2.^2);
        x = 1e4*[x1 x2];
        signaldegree = '2';
        signalregularity = 'approx. ---> 2-';
        wavelet = 'db3';
        maxlevel = 6;
        detectlevel = 1;
        
    case {7,'piecewise exponential signal with zero-order derivative discontinuity'}                   
    case {8,'piecewise exponential signal with first-order derivative discontinuity'}        
    case {9,'piecewise exponential signal with second-order derivative discontinuity'}     
        
    case {10,'two close singularities'}
        x = [400*(500:-1:1)/500 zeros(1,10) 400*(-1:-1:-514)/500 ];
        signaldegree = '1';
        signalregularity = 'approx. ---> 1-';
        wavelet = 'db2';
        maxlevel = 6;
        detectlevel = 1;
        
    case {11,'frequency discontinuity'}
        f1 = 87;
        f2 = 209;
        t = (0:N/4-1)*Ts;
        x1 = sin(2*pi*f1*t);
        x2 = sin(2*pi*f2*t);
        x = [x1 x2 x1 x2];
        signaldegree = 'Inf';
        signalregularity = '0';
        wavelet = 'db3';
        maxlevel = 6;
        detectlevel = 1;       
end

signallength = length(x);

% ----------------------------------------------------------------------------------------

% wavelet regularities and vanishing moments
switch wavelet
    case 'db1'
        waveletregularity = '0';
        vanishmoments = '1';
    case 'db2'
        waveletregularity = '0.50';
        vanishmoments = '2';
    case 'db3'
        waveletregularity = '0.91';
        vanishmoments = '3';
    case 'db4'
        waveletregularity = '1.27';
        vanishmoments = '4';
    case 'db5'
        waveletregularity = '1.59';
        vanishmoments = '5';
    case 'db6'
        waveletregularity = 'approx. 1.87';
        vanishmoments = '6';
    case 'db7'
        waveletregularity = '2.15';
        vanishmoments = '7';
    case 'db8'
        waveletregularity = 'approx. 2.40';
        vanishmoments = '8';
    case 'db9'
        waveletregularity = 'approx. 2.65';
        vanishmoments = '9';
    case 'db10'
        waveletregularity = '2.9';
        vanishmoments = '10';
end

% ---------------------------------------------------------------------------------------------

noisetype = 0;
switch noisetype
    case {0,'NONE'}
      noise = zeros(1,signallength);  
    case {1,'AWGN'}
      sigma = 0.2;
      noise = sigma*randn(1,signallength);  
    case {2,'UNIFORM'}
    case {3,'LPF-AWGN'}
    case {4,'NON-STATIONARY AWGN'}
end

% ---------------------------------------------------------------------------------------------

s = x + noise;

% ---------------------------------------------------------------------------------------------

% perform DWT/FWT (signal ---> coefficients cAJ,cDj)
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
    plot(cD(j,1:levellength(j)),'k-'),ylabel(['cD ' num2str(j)]),xlim([0 levellength(j)-1])
end
subplot(maxlevel+2,2^maxlevel,(maxlevel+1)*2^maxlevel+1)
plot(cA,'b-'),ylabel(['cA ' num2str(maxlevel)]),xlim([0 length(cA)-1])

% plot DFT/FWT coefficients expanded to signal length
figure('Name',['Signal and decomposition in DWT/FWT coefficients for analysis with Wavelet ' wavelet , '  (coefficients expanded to signal length)'])
subplot(maxlevel+2,1,1)
plot(s,'r-'),ylabel('signal'),axis([0 signallength-1 min(s) max(s)])
for j=1:maxlevel
    subplot(maxlevel+2,1,j+1)
    plot(cD(j,1:levellength(j)),'k-'),ylabel(['cD ' num2str(j)]),xlim([0 levellength(j)-1])
end
subplot(maxlevel+2,1,maxlevel+2)
plot(cA,'b-'),ylabel(['cA ' num2str(maxlevel)]),xlim([0 length(cA)-1])

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
plot(s,'r-'),ylabel('signal'),xlim([0 signallength-1]),title('approximations')
subplot(maxlevel+1,2,2)
plot(s,'r-'),ylabel('signal'),xlim([0 signallength-1]),title('details')
for j = 1:maxlevel
    subplot(maxlevel+1,2,2*j+1)
    plot(A(j,:),'b-'),ylabel(['A ' num2str(j)]),xlim([0 signallength-1])
    subplot(maxlevel+1,2,2*j+2)
    plot(D(j,:),'k-'),ylabel(['D ' num2str(j)]),xlim([0 signallength-1])
end

% ---------------------------------------------------------------------------------------------

% perform classical singularity detection for comparison
Nfilter = 100;
delay = Nfilter/2;
fc = 500;                       % signal/fc  1/1000, 2/500, 3/500, 4/1000, 5/1000, 6/500
h = fir1(Nfilter, fc/(fs/2),'high');
s_mod_filter = filter(h,1,s);
s_mod_filter = [s_mod_filter(delay+1:length(s)) zeros(1,delay)];
      
% ---------------------------------------------------------------------------------------------

% plot signal and detected singularities
figure('Name', 'Signal and detected singularities')     
subplot(3,1,1)
plot(n,s,'r-','LineWidth',3),title(['SIGNAL     (polynom.degree of signal = ' signaldegree ', regularity of signal = ' signalregularity ' )']),xlim([0 signallength-1])
subplot(3,1,2)
plot(n,D(detectlevel,:),'b-','LineWidth',2),title(['WAVELET analysis detected singularity    (wavelet = ' wavelet ', level = D' num2str(detectlevel) ', vanishing moments = ' vanishmoments ', regularity of wavelet = ' waveletregularity ')' ]),xlim([0 signallength-1]),ylim([1.2*min(D(detectlevel,50:950)) 1.2*max(D(detectlevel,50:950))])
subplot(3,1,3)
plot(n,s_mod_filter,'k-','LineWidth',2),title('CLASSICAL analysis (HPF) detected singularity'),xlim([0 signallength-1]),ylim([1.2*min(s_mod_filter(50:950)) 1.2*max(s_mod_filter(50:950))])

% ---------------------------------------------------------------------------------------------

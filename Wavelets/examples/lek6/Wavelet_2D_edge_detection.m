clc;clear;close all;
%-----------------------------------------------------------------------------------------------------
% BASIC EDGE DETECTION SCHEME FOR 2D WAVELET DSP:
% load/define image 
%    - grayscale images
%    - indexed (colormap) color pictures with smooth colormap
% analysis: perform DWT/FWT
% DSP: remove approximation coefficients
% synthesis: perform IDWT/IFWT on modified coefficients
% compare original and procesed images
%-----------------------------------------------------------------------------------------------------
% define or load signal                           
signaltype = 1;
switch signaltype
    case 1
        Iorig = imread('M:\MATLAB m-filer\testfigure.tif');
        Iorig = double(Iorig);
        Iorig(140:215,60:150) = 229*ones(76,91);
        [rows columns] = size(Iorig);
        nbcol = 256;
        colormap(gray(nbcol));
        map = colormap;
        close all;
end
%---------------------------------------------------------------------------------------------------------
% add noise
sigma = 20;
noise = sigma*randn(rows,columns);
I = Iorig + noise;
%---------------------------------------------------------------------------------------------------------
% perform DWT/FWT (signal ---> coefficients cAJ,cDj_h,cDj_v,cDj_d  1<=j<=J)
wavelet = 'sym4';
maxlevel = 1;
dwtmode('sym','nodisp');
[C,S] = wavedec2(I,maxlevel,wavelet);
%------------------------------------------------------------------------------------------------------------
% plot DFT/FWT coefficients at actual length
figure('Name',['ANALYSIS: signal and decomposition in DWT/FWT coefficients for analysis with wavelet = '...
    wavelet , ' at level = ' num2str(maxlevel) '    (coefficients shown at correct dyadic size ~ border distortion due to FWT removed)'])
subplot(1,2,1)
image(wcodemat(I,nbcol));colormap(map);colorbar;axis image,title(['signal, size = [' num2str(rows) ' ' num2str(columns) ']'])
subplot(1,2,2)
j = 1:maxlevel;
lines_x_coord = [columns./2.^j+1 zeros(1,maxlevel); columns./2.^j+1 columns./2.^(j-1)];
lines_y_coord = [zeros(1,maxlevel) rows./2.^j+1; rows./2.^(j-1) rows./2.^j+1];
for j = maxlevel:-1:1
    if j == maxlevel
        cA = appcoef2(C,S,wavelet,maxlevel);
        row_shrinkage = S(1,1) - rows/2^maxlevel;
        row_shrinkage_top = floor(row_shrinkage/2);
        column_shrinkage = S(1,2) - columns/2^maxlevel;
        column_shrinkage_left = floor(column_shrinkage/2);
        cA_central = cA((row_shrinkage_top + 1):(row_shrinkage_top + rows/2^maxlevel),(column_shrinkage_left + 1):(column_shrinkage_left + columns/2^maxlevel)); 
    end
    cDh = detcoef2('h',C,S,j);
    row_shrinkage = S(2 + maxlevel - j,1) - rows/2^j;
    row_shrinkage_top = floor(row_shrinkage/2);
    column_shrinkage = S(2 + maxlevel - j,2) - columns/2^j;
    column_shrinkage_left = floor(column_shrinkage/2);
    cDh_central = cDh((row_shrinkage_top + 1):(row_shrinkage_top + rows/2^j),(column_shrinkage_left + 1):(column_shrinkage_left + columns/2^j)); 
    cDv = detcoef2('v',C,S,j);
    cDv_central = cDv((row_shrinkage_top + 1):(row_shrinkage_top + rows/2^j),(column_shrinkage_left + 1):(column_shrinkage_left + columns/2^j)); 
    cDd = detcoef2('d',C,S,j);
    cDd_central = cDd((row_shrinkage_top + 1):(row_shrinkage_top + rows/2^j),(column_shrinkage_left + 1):(column_shrinkage_left + columns/2^j)); 
    if j == maxlevel
        Icoeffs = [wcodemat(cA_central,nbcol) wcodemat(cDh_central,nbcol); wcodemat(cDv_central,nbcol) wcodemat(cDd_central,nbcol)];
    else
        Icoeffs = [Icoeffs wcodemat(cDh_central,nbcol); wcodemat(cDv_central,nbcol) wcodemat(cDd_central,nbcol)];
    end
end
image(wcodemat(Icoeffs,nbcol));colormap(map);colorbar;axis image,title(['cA: approx. coeffs [' num2str(rows) '/2^l^e^v^e^l ' num2str(columns) '/2^l^e^v^e^l]                  cDh: horisontal detail coeffs [' num2str(rows) '/2^l^e^v^e^l ' num2str(columns) '/2^l^e^v^e^l]']),...
    xlabel(['cDv: vertical detail coeffs [' num2str(rows) '/2^l^e^v^e^l ' num2str(columns) '/2^l^e^v^e^l]           cDd: diagonal detail coeffs [' num2str(rows) '/2^l^e^v^e^l ' num2str(columns) '/2^l^e^v^e^l]']),line(lines_x_coord,lines_y_coord,'LineWidth',2,'Color','b')
%------------------------------------------------------------------------------------------------------------
% plot DFT/FWT coefficients expanded to signal length
figure('Name',['ANALYSIS: signal and decomposition in DWT/FWT coefficients for analysis with wavelet = '...
    wavelet , '    (complete coefficients including border distortion due to FWT, shown at signal size)'])
subplot(maxlevel+1,4,1)
image(wcodemat(I,nbcol));colormap(map);colorbar;axis image,title(['signal, size = [' num2str(rows) ' ' num2str(columns) ']'])
for j = 1:maxlevel
    subplot(maxlevel+1,4,j*4+1)
    cA = appcoef2(C,S,wavelet,j);
    image(wcodemat(cA,nbcol));colormap(map);colorbar;axis image,ylabel(['level = ' num2str(j)']),title(['cA: approx. coeffs [' num2str(size(cA,1)) ' ' num2str(size(cA,2)) ']'])
    subplot(maxlevel+1,4,j*4+2)
    cDh = detcoef2('h',C,S,j);
    image(wcodemat(cDh,nbcol));colormap(map);colorbar;axis image,title(['cDh: horisontal detail coeffs [' num2str(size(cDh,1)) ' ' num2str(size(cDh,2)) ']'])
    subplot(maxlevel+1,4,j*4+3)
    cDv = detcoef2('v',C,S,j);
    image(wcodemat(cDv,nbcol));colormap(map);colorbar;axis image,title(['cDv: vertical detail coeffs [' num2str(size(cDv,1)) ' ' num2str(size(cDv,2)) ']'])
    subplot(maxlevel+1,4,j*4+4)
    cDd = detcoef2('d',C,S,j);
    image(wcodemat(cDd,nbcol));colormap(map);colorbar;axis image,title(['cDd: diagonal detail coeffs [' num2str(size(cDd,1)) ' ' num2str(size(cDd,2)) ']'])
end
%-----------------------------------------------------------------------------------------------------------
% plot reconstructed approximations and details
figure('Name',['SYNTHESIS: signal and contributes from individually approximations and details reconstructed with wavelet = '...
    wavelet ' and truncated to correct signal size ~ border distortion due to IFWT removed'])
subplot(maxlevel+1,4,1)
image(wcodemat(I,nbcol));colormap(map);colorbar;axis image,title(['signal, size = [' num2str(rows) ' ' num2str(columns) ']'])
for j = 1:maxlevel
    subplot(maxlevel+1,4,j*4+1)
    A = wrcoef2('a',C,S,wavelet,j);
    image(wcodemat(A,nbcol));colormap(map);colorbar;axis image,ylabel(['level = ' num2str(j)']),title(['A: approximation [' num2str(rows) ' ' num2str(columns) ']'])
    subplot(maxlevel+1,4,j*4+2)
    Dh = wrcoef2('h',C,S,wavelet,j);
    image(wcodemat(Dh,nbcol));colormap(map);colorbar;axis image,title(['Dh: horisontal detail [' num2str(rows) ' ' num2str(columns) ']'])
    subplot(maxlevel+1,4,j*4+3)
    Dv = wrcoef2('v',C,S,wavelet,j);
    image(wcodemat(Dv,nbcol));colormap(map);colorbar;axis image,title(['Dv: vertical detail [' num2str(rows) ' ' num2str(columns) ']'])
    subplot(maxlevel+1,4,j*4+4)
    Dd = wrcoef2('d',C,S,wavelet,j);
    image(wcodemat(Dd,nbcol));colormap(map);colorbar;axis image,title(['Dd: diagonal detail [' num2str(rows) ' ' num2str(columns) ']'])
end
%-------------------------------------------------------------------------------------------------------------
% DSP: modify coefficients
Cmod = C;
edge_detection_method = 1;
switch edge_detection_method
    case {0,'NONE'}
                      
    case {1,'REMOVE APPROXIMATION COEFFICIENTS AT DEEPEST LEVEL'}
        cA = appcoef2(Cmod,S,wavelet,maxlevel);
        cA(:,:) = 0;
        Cmod = modify_2D_wavelet_coeffs(Cmod,S,'a',maxlevel,cA);
end
Imod = waverec2(Cmod,S,wavelet);
%------------------------------------------------------------------------------------------------------------
% plot modified DFT/FWT coefficients at actual length
figure('Name',['ANALYSIS: processed signal and decomposition in DWT/FWT coefficients for analysis with wavelet = '...
    wavelet , ' at level = ' num2str(maxlevel) '    (coefficients shown at correct dyadic size ~ border distortion due to FWT removed)'])
subplot(1,2,1)
image(wcodemat(Imod,nbcol));colormap(map);colorbar;axis image,title(['signal, size = [' num2str(rows) ' ' num2str(columns) ']'])
subplot(1,2,2)
j = 1:maxlevel;
lines_x_coord = [columns./2.^j+1 zeros(1,maxlevel); columns./2.^j+1 columns./2.^(j-1)];
lines_y_coord = [zeros(1,maxlevel) rows./2.^j+1; rows./2.^(j-1) rows./2.^j+1];
for j = maxlevel:-1:1
    if j == maxlevel
        cA = appcoef2(Cmod,S,wavelet,maxlevel);
        row_shrinkage = S(1,1) - rows/2^maxlevel;
        row_shrinkage_top = floor(row_shrinkage/2);
        column_shrinkage = S(1,2) - columns/2^maxlevel;
        column_shrinkage_left = floor(column_shrinkage/2);
        cA_central = cA((row_shrinkage_top + 1):(row_shrinkage_top + rows/2^maxlevel),(column_shrinkage_left + 1):(column_shrinkage_left + columns/2^maxlevel)); 
    end
    cDh = detcoef2('h',Cmod,S,j);
    row_shrinkage = S(2 + maxlevel - j,1) - rows/2^j;
    row_shrinkage_top = floor(row_shrinkage/2);
    column_shrinkage = S(2 + maxlevel - j,2) - columns/2^j;
    column_shrinkage_left = floor(column_shrinkage/2);
    cDh_central = cDh((row_shrinkage_top + 1):(row_shrinkage_top + rows/2^j),(column_shrinkage_left + 1):(column_shrinkage_left + columns/2^j)); 
    cDv = detcoef2('v',Cmod,S,j);
    cDv_central = cDv((row_shrinkage_top + 1):(row_shrinkage_top + rows/2^j),(column_shrinkage_left + 1):(column_shrinkage_left + columns/2^j)); 
    cDd = detcoef2('d',Cmod,S,j);
    cDd_central = cDd((row_shrinkage_top + 1):(row_shrinkage_top + rows/2^j),(column_shrinkage_left + 1):(column_shrinkage_left + columns/2^j)); 
    if j == maxlevel
        Icoeffs = [wcodemat(cA_central,nbcol) wcodemat(cDh_central,nbcol); wcodemat(cDv_central,nbcol) wcodemat(cDd_central,nbcol)];
    else
        Icoeffs = [Icoeffs wcodemat(cDh_central,nbcol); wcodemat(cDv_central,nbcol) wcodemat(cDd_central,nbcol)];
    end
end
image(wcodemat(Icoeffs,nbcol));colormap(map);colorbar;axis image,title(['cA: approx. coeffs [' num2str(rows) '/2^l^e^v^e^l ' num2str(columns) '/2^l^e^v^e^l]                  cDh: horisontal detail coeffs [' num2str(rows) '/2^l^e^v^e^l ' num2str(columns) '/2^l^e^v^e^l]']),...
    xlabel(['cDv: vertical detail coeffs [' num2str(rows) '/2^l^e^v^e^l ' num2str(columns) '/2^l^e^v^e^l]           cDd: diagonal detail coeffs [' num2str(rows) '/2^l^e^v^e^l ' num2str(columns) '/2^l^e^v^e^l]']),line(lines_x_coord,lines_y_coord,'LineWidth',2,'Color','b')
%------------------------------------------------------------------------------------------------------------
% plot modified DFT/FWT coefficients expanded to signal length
figure('Name',['ANALYSIS: processed signal and decomposition in DWT/FWT coefficients for analysis with wavelet = '...
    wavelet , '    (complete coefficients including border distortion due to FWT, shown at signal size)'])
subplot(maxlevel+1,4,1)
image(wcodemat(Imod,nbcol));colormap(map);colorbar;axis image,title(['signal, size = [' num2str(rows) ' ' num2str(columns) ']'])
for j = 1:maxlevel
    subplot(maxlevel+1,4,j*4+1)
    cA = appcoef2(Cmod,S,wavelet,j);
    image(wcodemat(cA,nbcol));colormap(map);colorbar;axis image,ylabel(['level = ' num2str(j)']),title(['cA: approx. coeffs [' num2str(size(cA,1)) ' ' num2str(size(cA,2)) ']'])
    subplot(maxlevel+1,4,j*4+2)
    cDh = detcoef2('h',Cmod,S,j);
    image(wcodemat(cDh,nbcol));colormap(map);colorbar;axis image,title(['cDh: horisontal detail coeffs [' num2str(size(cDh,1)) ' ' num2str(size(cDh,2)) ']'])
    subplot(maxlevel+1,4,j*4+3)
    cDv = detcoef2('v',Cmod,S,j);
    image(wcodemat(cDv,nbcol));colormap(map);colorbar;axis image,title(['cDv: vertical detail coeffs [' num2str(size(cDv,1)) ' ' num2str(size(cDv,2)) ']'])
    subplot(maxlevel+1,4,j*4+4)
    cDd = detcoef2('d',Cmod,S,j);
    image(wcodemat(cDd,nbcol));colormap(map);colorbar;axis image,title(['cDd: diagonal detail coeffs [' num2str(size(cDd,1)) ' ' num2str(size(cDd,2)) ']'])
end
%-----------------------------------------------------------------------------------------------------------
% plot modified reconstructed approximations and details
figure('Name',['SYNTHESIS: processed signal and contributes from individually approximations and details reconstructed with wavelet = '...
    wavelet ' and truncated to correct signal size ~ border distortion due to IFWT removed'])
subplot(maxlevel+1,4,1)
image(wcodemat(Imod,nbcol));colormap(map);colorbar;axis image,title(['signal, size = [' num2str(rows) ' ' num2str(columns) ']'])
for j = 1:maxlevel
    subplot(maxlevel+1,4,j*4+1)
    A = wrcoef2('a',Cmod,S,wavelet,j);
    image(wcodemat(A,nbcol));colormap(map);colorbar;axis image,ylabel(['level = ' num2str(j)']),title(['A: approximation [' num2str(rows) ' ' num2str(columns) ']'])
    subplot(maxlevel+1,4,j*4+2)
    Dh = wrcoef2('h',Cmod,S,wavelet,j);
    image(wcodemat(Dh,nbcol));colormap(map);colorbar;axis image,title(['Dh: horisontal detail [' num2str(rows) ' ' num2str(columns) ']'])
    subplot(maxlevel+1,4,j*4+3)
    Dv = wrcoef2('v',Cmod,S,wavelet,j);
    image(wcodemat(Dv,nbcol));colormap(map);colorbar;axis image,title(['Dv: vertical detail [' num2str(rows) ' ' num2str(columns) ']'])
    subplot(maxlevel+1,4,j*4+4)
    Dd = wrcoef2('d',Cmod,S,wavelet,j);
    image(wcodemat(Dd,nbcol));colormap(map);colorbar;axis image,title(['Dd: diagonal detail [' num2str(rows) ' ' num2str(columns) ']'])
end
% ---------------------------------------------------------------------------------------------
%plot image and processed image
figure('Name',['Original image, noisy image and processed image'])
subplot(1,3,1)
image(wcodemat(Iorig,nbcol));colormap(map);colorbar;axis image,title(['original signal, size = [' num2str(rows) ' ' num2str(columns) ']'])
subplot(1,3,2)
image(wcodemat(I,nbcol));colormap(map);colorbar;axis image,title(['noisy signal, size = [' num2str(rows) ' ' num2str(columns) ']']),...
subplot(1,3,3)
edges = find(abs(Imod) >= 65);
Imod(edges) = 255;
notedges = find(abs(Imod) < 65);
Imod(notedges) = 0;
image(wcodemat(Imod,nbcol));colormap(map);colorbar;axis image,title(['processed signal, size = [' num2str(rows) ' ' num2str(columns) ']']),...
figure('Name',['Original image'])
image(wcodemat(Iorig,nbcol));colormap(map);colorbar;axis image
figure('Name',['Noisy image'])
image(wcodemat(I,nbcol));colormap(map);colorbar;axis image
figure('Name',['Wavelet processed image        NB: NOT binary output'])
image(wcodemat(Imod,nbcol));colormap(map);colorbar;axis image
%--------------------------------------------------------------------------------------------------------------
% comparison with standard edge detection schemes
figure('Name','Classic edge detection with Sobel algorithm (thresholded gradient based)     NB: binary output')
Isobel=edge(I,'sobel');
imshow(Isobel);colormap(map);colorbar;axis image
figure('Name','Classic edge detection with Canny algorithm (second derivative zerocrossings based)    NB: binary output')
Icanny=edge(I,'canny');
imshow(Icanny);colormap(map);colorbar;axis image
%---------------------------------------------------------------------------------------------------------------
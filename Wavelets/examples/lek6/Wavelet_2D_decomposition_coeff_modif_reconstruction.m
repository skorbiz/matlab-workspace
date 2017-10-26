clc;clear;close all;
%-----------------------------------------------------------------------------------------------------
% BASIC SCHEME FOR 2D WAVELET DSP:
% load/define image 
%    - grayscale images
%    - indexed (colormap) color pictures with smooth colormap
% analysis: perform DWT/FWT
% DSP: modify coefficients
% synthesis: perform IDWT/IFWT on modified coefficients
% compare original and processed images
%-----------------------------------------------------------------------------------------------------
% define or load signal                           
signaltype = 1;
switch signaltype
    case 1
        load wbarb;  % --> X = wbarb; map; indexed image
        Iorig = X;
        [rows columns] = size(Iorig);
        nbcol = size(map,1);
        
    case 2
        load detfingr;  % --> X = detfingr; map; indexed image
        Iorig = X;
        [rows columns] = size(Iorig);
        nbcol = 192;
        colormap(gray(nbcol));
        map = colormap; 
        close all;
                
    case 3
        load detail;  % --> X = detail; map; indexed image
        Iorig = X(5:356,2:369);
        [rows columns] = size(Iorig);
        nbcol = 64;
        colormap(pink(nbcol));
        map = colormap;
        close all;
end
%---------------------------------------------------------------------------------------------------------
% add noise
sigma = 45;
noise = sigma*randn(rows,columns);
I = Iorig + noise;
%------------------------------------------------------------------------------------------------------------
% perform DWT/FWT (signal ---> coefficients cAJ,cDj_h,cDj_v,cDj_d  1<=j<=J)
wavelet = 'bior3.7';
maxlevel = 3;
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
SP_type = 1;
switch SP_type
    case {1, 'treat h,v,d equally at each level'}
        % upper left part of image: do nothing
        % upper right part of image: remove cD1 (h,v,d)
        cDh = detcoef2('h',Cmod,S,1);
        [coeffrows1 coeffcolumns1] = size(cDh);
        cDh(1:floor(coeffrows1/2),floor(coeffcolumns1/2):coeffcolumns1) = zeros(floor(coeffrows1/2),coeffcolumns1 - floor(coeffcolumns1/2) + 1);
        Cmod = modify_2D_wavelet_coeffs(Cmod,S,'h',1,cDh);
        cDv = detcoef2('v',Cmod,S,1);
        cDv(1:floor(coeffrows1/2),floor(coeffcolumns1/2):coeffcolumns1) = zeros(floor(coeffrows1/2),coeffcolumns1 - floor(coeffcolumns1/2) + 1);
        Cmod = modify_2D_wavelet_coeffs(Cmod,S,'v',1,cDv);
        cDd = detcoef2('d',Cmod,S,1);
        cDd(1:floor(coeffrows1/2),floor(coeffcolumns1/2):coeffcolumns1) = zeros(floor(coeffrows1/2),coeffcolumns1 - floor(coeffcolumns1/2) + 1);
        Cmod = modify_2D_wavelet_coeffs(Cmod,S,'d',1,cDd);
        % lower left part of image: remove cD1, cD2 (h,v,d)
        cDh = detcoef2('h',Cmod,S,1);
        [coeffrows1 coeffcolumns1] = size(cDh);
        cDh(floor(coeffrows1/2):coeffrows1,1:floor(coeffcolumns1/2)) = zeros(coeffrows1 - floor(coeffrows1/2) + 1, floor(coeffcolumns1/2));
        Cmod = modify_2D_wavelet_coeffs(Cmod,S,'h',1,cDh);
        cDv = detcoef2('v',Cmod,S,1);
        cDv(floor(coeffrows1/2):coeffrows1,1:floor(coeffcolumns1/2)) = zeros(coeffrows1 - floor(coeffrows1/2) + 1, floor(coeffcolumns1/2));
        Cmod = modify_2D_wavelet_coeffs(Cmod,S,'v',1,cDv);
        cDd = detcoef2('d',Cmod,S,1);
        cDd(floor(coeffrows1/2):coeffrows1,1:floor(coeffcolumns1/2)) = zeros(coeffrows1 - floor(coeffrows1/2) + 1, floor(coeffcolumns1/2));
        Cmod = modify_2D_wavelet_coeffs(Cmod,S,'d',1,cDd);
        cDh = detcoef2('h',Cmod,S,2);
        [coeffrows2 coeffcolumns2] = size(cDh);
        cDh(floor(coeffrows2/2):coeffrows2,1:floor(coeffcolumns2/2)) = zeros(coeffrows2 - floor(coeffrows2/2) + 1, floor(coeffcolumns2/2));
        Cmod = modify_2D_wavelet_coeffs(Cmod,S,'h',2,cDh);
        cDv = detcoef2('v',Cmod,S,2);
        cDv(floor(coeffrows2/2):coeffrows2,1:floor(coeffcolumns2/2)) = zeros(coeffrows2 - floor(coeffrows2/2) + 1, floor(coeffcolumns2/2));
        Cmod = modify_2D_wavelet_coeffs(Cmod,S,'v',2,cDv);
        cDd = detcoef2('d',Cmod,S,2);
        cDd(floor(coeffrows2/2):coeffrows2,1:floor(coeffcolumns2/2)) = zeros(coeffrows2 - floor(coeffrows2/2) + 1, floor(coeffcolumns2/2));
        Cmod = modify_2D_wavelet_coeffs(Cmod,S,'d',2,cDd);
        % lower right part of image: remove cD1, cD2, cD3 (h,v,d)
        cDh = detcoef2('h',Cmod,S,1);
        [coeffrows1 coeffcolumns1] = size(cDh);
        cDh(floor(coeffrows1/2):coeffrows1,floor(coeffcolumns1/2):coeffcolumns1) = zeros(coeffrows1 - floor(coeffrows1/2) + 1, coeffcolumns1 - floor(coeffcolumns1/2) + 1);
        Cmod = modify_2D_wavelet_coeffs(Cmod,S,'h',1,cDh);
        cDv = detcoef2('v',Cmod,S,1);
        cDv(floor(coeffrows1/2):coeffrows1,floor(coeffcolumns1/2):coeffcolumns1) = zeros(coeffrows1 - floor(coeffrows1/2) + 1, coeffcolumns1 - floor(coeffcolumns1/2) + 1);
        Cmod = modify_2D_wavelet_coeffs(Cmod,S,'v',1,cDv);
        cDd = detcoef2('d',Cmod,S,1);
        cDd(floor(coeffrows1/2):coeffrows1,floor(coeffcolumns1/2):coeffcolumns1) = zeros(coeffrows1 - floor(coeffrows1/2) + 1, coeffcolumns1 - floor(coeffcolumns1/2) + 1);
        Cmod = modify_2D_wavelet_coeffs(Cmod,S,'d',1,cDd);
        cDh = detcoef2('h',Cmod,S,2);
        [coeffrows2 coeffcolumns2] = size(cDh);
        cDh(floor(coeffrows2/2):coeffrows2,floor(coeffcolumns2/2):coeffcolumns2) = zeros(coeffrows2 - floor(coeffrows2/2) + 1, coeffcolumns2 - floor(coeffcolumns2/2) + 1);
        Cmod = modify_2D_wavelet_coeffs(Cmod,S,'h',2,cDh);
        cDv = detcoef2('v',Cmod,S,2);
        cDv(floor(coeffrows2/2):coeffrows2,floor(coeffcolumns2/2):coeffcolumns2) = zeros(coeffrows2 - floor(coeffrows2/2) + 1, coeffcolumns2 - floor(coeffcolumns2/2) + 1);
        Cmod = modify_2D_wavelet_coeffs(Cmod,S,'v',2,cDv);
        cDd = detcoef2('d',Cmod,S,2);
        cDd(floor(coeffrows2/2):coeffrows2,floor(coeffcolumns2/2):coeffcolumns2) = zeros(coeffrows2 - floor(coeffrows2/2) + 1, coeffcolumns2 - floor(coeffcolumns2/2) + 1);
        Cmod = modify_2D_wavelet_coeffs(Cmod,S,'d',2,cDd);
        cDh = detcoef2('h',Cmod,S,3);
        [coeffrows3 coeffcolumns3] = size(cDh);
        cDh(floor(coeffrows3/2):coeffrows3,floor(coeffcolumns3/2):coeffcolumns3) = zeros(coeffrows3 - floor(coeffrows3/2) + 1, coeffcolumns3 - floor(coeffcolumns3/2) + 1);
        Cmod = modify_2D_wavelet_coeffs(Cmod,S,'h',3,cDh);
        cDv = detcoef2('v',Cmod,S,3);
        cDv(floor(coeffrows3/2):coeffrows3,floor(coeffcolumns3/2):coeffcolumns3) = zeros(coeffrows3 - floor(coeffrows3/2) + 1, coeffcolumns3 - floor(coeffcolumns3/2) + 1);
        Cmod = modify_2D_wavelet_coeffs(Cmod,S,'v',3,cDv);
        cDd = detcoef2('d',Cmod,S,3);
        cDd(floor(coeffrows3/2):coeffrows3,floor(coeffcolumns3/2):coeffcolumns3) = zeros(coeffrows3 - floor(coeffrows3/2) + 1, coeffcolumns3 - floor(coeffcolumns3/2) + 1);
        Cmod = modify_2D_wavelet_coeffs(Cmod,S,'d',3,cDd);
         
    case {2, 'treat h,v,d differently'}
        % upper left part of image: do nothing
        % upper right part of image: remove cD1h, cD2h 
        cDh = detcoef2('h',Cmod,S,1);
        [coeffrows1 coeffcolumns1] = size(cDh);
        cDh(1:floor(coeffrows1/2),floor(coeffcolumns1/2):coeffcolumns1) = zeros(floor(coeffrows1/2),coeffcolumns1 - floor(coeffcolumns1/2) + 1);
        Cmod = modify_2D_wavelet_coeffs(Cmod,S,'h',1,cDh);
        cDh = detcoef2('h',Cmod,S,2);
        [coeffrows2 coeffcolumns2] = size(cDh);
        cDh(1:floor(coeffrows2/2),floor(coeffcolumns2/2):coeffcolumns2) = zeros(floor(coeffrows2/2),coeffcolumns2 - floor(coeffcolumns2/2) + 1);
        Cmod = modify_2D_wavelet_coeffs(Cmod,S,'h',2,cDh);
        % lower left part of image: remove cD1v, cD2v
        cDv = detcoef2('v',Cmod,S,1);
        [coeffrows1 coeffcolumns1] = size(cDv);
        cDv(floor(coeffrows1/2):coeffrows1,1:floor(coeffcolumns1/2)) = zeros(coeffrows1 - floor(coeffrows1/2) + 1, floor(coeffcolumns1/2));
        Cmod = modify_2D_wavelet_coeffs(Cmod,S,'v',1,cDv);
        cDv = detcoef2('v',Cmod,S,2);
        [coeffrows2 coeffcolumns2] = size(cDv);
        cDv(floor(coeffrows2/2):coeffrows2,1:floor(coeffcolumns2/2)) = zeros(coeffrows2 - floor(coeffrows2/2) + 1, floor(coeffcolumns2/2));
        Cmod = modify_2D_wavelet_coeffs(Cmod,S,'v',2,cDv);
        % lower right part of image: remove cD1d, cD2d
        cDd = detcoef2('d',Cmod,S,1);
        [coeffrows1 coeffcolumns1] = size(cDd);
        cDd(floor(coeffrows1/2):coeffrows1,floor(coeffcolumns1/2):coeffcolumns1) = zeros(coeffrows1 - floor(coeffrows1/2) + 1, coeffcolumns1 - floor(coeffcolumns1/2) + 1);
        Cmod = modify_2D_wavelet_coeffs(Cmod,S,'d',1,cDd);
        cDd = detcoef2('d',Cmod,S,2);
        [coeffrows2 coeffcolumns2] = size(cDd);
        cDd(floor(coeffrows2/2):coeffrows2,floor(coeffcolumns2/2):coeffcolumns2) = zeros(coeffrows2 - floor(coeffrows2/2) + 1, coeffcolumns2 - floor(coeffcolumns2/2) + 1);
        Cmod = modify_2D_wavelet_coeffs(Cmod,S,'d',2,cDd);
end
Imod = waverec2(Cmod,S,wavelet);
%------------------------------------------------------------------------------------------------------------
% plot modified DFT/FWT coefficients at actual length
figure('Name',['ANALYSIS: modified signal and decomposition in DWT/FWT coefficients for analysis with wavelet = '...
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
figure('Name',['ANALYSIS: modified signal and decomposition in DWT/FWT coefficients for analysis with wavelet = '...
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
figure('Name',['SYNTHESIS: modified signal and contributes from individually approximations and details reconstructed with wavelet = '...
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
%--------------------------------------------------------------------------------------------------------------
%plot image and processed image
figure('Name',['Original image, noisy image and modified image'])
subplot(1,3,1)
image(wcodemat(Iorig,nbcol));colormap(map);colorbar;axis image,title(['original signal, size = [' num2str(rows) ' ' num2str(columns) ']'])
subplot(1,3,2)
image(wcodemat(I,nbcol));colormap(map);colorbar;axis image,title(['noisy signal, size = [' num2str(rows) ' ' num2str(columns) ']'])
subplot(1,3,3)
image(wcodemat(Imod,nbcol));colormap(map);colorbar;axis image,title(['modified signal, size = [' num2str(rows) ' ' num2str(columns) ']'])
%--------------------------------------------------------------------------------------------------------------

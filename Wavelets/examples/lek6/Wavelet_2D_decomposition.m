clc;clear;close all;
%-------------------------------------------------------------------------
% Wavelet decomposition of
% - grayscale images
% - indexed (colormap) color pictures with smooth colormap

% For truecolor RGB pictures (no colormap applied) and
% indexed pictures with non-smooth colormap,
% convert to grayscale before FWT
% (se special m-file)
%-------------------------------------------------------------------------
% define or load signal                           
signaltype = 1;
switch signaltype
    case 1
        load wbarb;  % --> X = wbarb; map; indexed image
        I = X;
        [rows columns] = size(I);
        nbcol = size(map,1);
        
    case 2
        load detfingr;  % --> X = detfingr; map; indexed image
        I = X;
        [rows columns] = size(I);
        nbcol = 192;
        colormap(gray(nbcol));
        map = colormap; 
        close all;
                
    case 3
        load detail;  % --> X = detail; map; indexed image
        I = X(5:356,2:369);
        [rows columns] = size(I);
        nbcol = 64;
        colormap(pink(nbcol));
        map = colormap;
        close all;
end
%------------------------------------------------------------------------------------------------------------
% perform DWT/FWT (signal ---> coefficients cAJ,cDj_h,cDj_v,cDj_d  1<=j<=J)
wavelet = 'haar';
maxlevel = 2;
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

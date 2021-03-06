clc;clear;close all;
%-----------------------------------------------------------------------------------------------------
% BASIC COMPRESSION SCHEMES FOR 2D WAVELET DSP:
% load/define image 
%    - grayscale images
%    - indexed (colormap) color pictures with smooth colormap
% analysis: perform DWT/FWT
% DSP: perform global or level-wise thresholding of coefficients
% synthesis: perform IDWT/IFWT on modified coefficients
% compare original and compressed images
%-----------------------------------------------------------------------------------------------------
% define or load signal                           
signaltype = 2;
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
%---------------------------------------------------------------------------------------------------------
% perform DWT/FWT (signal ---> coefficients cAJ,cDj_h,cDj_v,cDj_d  1<=j<=J)
wavelet = 'bior3.7';
compresslevel = 3;
dwtmode('sym','nodisp');
[C,S] = wavedec2(I,compresslevel,wavelet);
%------------------------------------------------------------------------------------------------------------
% plot DFT/FWT coefficients at actual length
figure('Name',['ANALYSIS: signal and decomposition in DWT/FWT coefficients for analysis with wavelet = '...
    wavelet , ' at level = ' num2str(compresslevel) '    (coefficients shown at correct dyadic size ~ border distortion due to FWT removed)'])
subplot(1,2,1)
image(wcodemat(I,nbcol));colormap(map);colorbar;axis image,title(['signal, size = [' num2str(rows) ' ' num2str(columns) ']'])
subplot(1,2,2)
j = 1:compresslevel;
lines_x_coord = [columns./2.^j+1 zeros(1,compresslevel); columns./2.^j+1 columns./2.^(j-1)];
lines_y_coord = [zeros(1,compresslevel) rows./2.^j+1; rows./2.^(j-1) rows./2.^j+1];
for j = compresslevel:-1:1
    if j == compresslevel
        cA = appcoef2(C,S,wavelet,compresslevel);
        row_shrinkage = S(1,1) - rows/2^compresslevel;
        row_shrinkage_top = floor(row_shrinkage/2);
        column_shrinkage = S(1,2) - columns/2^compresslevel;
        column_shrinkage_left = floor(column_shrinkage/2);
        cA_central = cA((row_shrinkage_top + 1):(row_shrinkage_top + rows/2^compresslevel),(column_shrinkage_left + 1):(column_shrinkage_left + columns/2^compresslevel)); 
    end
    cDh = detcoef2('h',C,S,j);
    row_shrinkage = S(2 + compresslevel - j,1) - rows/2^j;
    row_shrinkage_top = floor(row_shrinkage/2);
    column_shrinkage = S(2 + compresslevel - j,2) - columns/2^j;
    column_shrinkage_left = floor(column_shrinkage/2);
    cDh_central = cDh((row_shrinkage_top + 1):(row_shrinkage_top + rows/2^j),(column_shrinkage_left + 1):(column_shrinkage_left + columns/2^j)); 
    cDv = detcoef2('v',C,S,j);
    cDv_central = cDv((row_shrinkage_top + 1):(row_shrinkage_top + rows/2^j),(column_shrinkage_left + 1):(column_shrinkage_left + columns/2^j)); 
    cDd = detcoef2('d',C,S,j);
    cDd_central = cDd((row_shrinkage_top + 1):(row_shrinkage_top + rows/2^j),(column_shrinkage_left + 1):(column_shrinkage_left + columns/2^j)); 
    if j == compresslevel
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
subplot(compresslevel+1,4,1)
image(wcodemat(I,nbcol));colormap(map);colorbar;axis image,title(['signal, size = [' num2str(rows) ' ' num2str(columns) ']'])
for j = 1:compresslevel
    subplot(compresslevel+1,4,j*4+1)
    cA = appcoef2(C,S,wavelet,j);
    image(wcodemat(cA,nbcol));colormap(map);colorbar;axis image,ylabel(['level = ' num2str(j)']),title(['cA: approx. coeffs [' num2str(size(cA,1)) ' ' num2str(size(cA,2)) ']'])
    subplot(compresslevel+1,4,j*4+2)
    cDh = detcoef2('h',C,S,j);
    image(wcodemat(cDh,nbcol));colormap(map);colorbar;axis image,title(['cDh: horisontal detail coeffs [' num2str(size(cDh,1)) ' ' num2str(size(cDh,2)) ']'])
    subplot(compresslevel+1,4,j*4+3)
    cDv = detcoef2('v',C,S,j);
    image(wcodemat(cDv,nbcol));colormap(map);colorbar;axis image,title(['cDv: vertical detail coeffs [' num2str(size(cDv,1)) ' ' num2str(size(cDv,2)) ']'])
    subplot(compresslevel+1,4,j*4+4)
    cDd = detcoef2('d',C,S,j);
    image(wcodemat(cDd,nbcol));colormap(map);colorbar;axis image,title(['cDd: diagonal detail coeffs [' num2str(size(cDd,1)) ' ' num2str(size(cDd,2)) ']'])
end
%-----------------------------------------------------------------------------------------------------------
% plot reconstructed approximations and details
figure('Name',['SYNTHESIS: signal and contributes from individually approximations and details reconstructed with wavelet = '...
    wavelet ' and truncated to correct signal size ~ border distortion due to IFWT removed'])
subplot(compresslevel+1,4,1)
image(wcodemat(I,nbcol));colormap(map);colorbar;axis image,title(['signal, size = [' num2str(rows) ' ' num2str(columns) ']'])
for j = 1:compresslevel
    subplot(compresslevel+1,4,j*4+1)
    A = wrcoef2('a',C,S,wavelet,j);
    image(wcodemat(A,nbcol));colormap(map);colorbar;axis image,ylabel(['level = ' num2str(j)']),title(['A: approximation [' num2str(rows) ' ' num2str(columns) ']'])
    subplot(compresslevel+1,4,j*4+2)
    Dh = wrcoef2('h',C,S,wavelet,j);
    image(wcodemat(Dh,nbcol));colormap(map);colorbar;axis image,title(['Dh: horisontal detail [' num2str(rows) ' ' num2str(columns) ']'])
    subplot(compresslevel+1,4,j*4+3)
    Dv = wrcoef2('v',C,S,wavelet,j);
    image(wcodemat(Dv,nbcol));colormap(map);colorbar;axis image,title(['Dv: vertical detail [' num2str(rows) ' ' num2str(columns) ']'])
    subplot(compresslevel+1,4,j*4+4)
    Dd = wrcoef2('d',C,S,wavelet,j);
    image(wcodemat(Dd,nbcol));colormap(map);colorbar;axis image,title(['Dd: diagonal detail [' num2str(rows) ' ' num2str(columns) ']'])
end
%-------------------------------------------------------------------------------------------------------------
% DSP: modify coefficients
Cmod = C;
compression_method = 2;
switch compression_method
    case {0,'NONE'}
            
    case {1,'GLOBAL HARD THRESHOLDING OF H,V,D DETAIL COEFFS'}
        switch signaltype
            case 1                      % supply one global threshold for level 1 to compresslevel
                threshold = 135.4;  %level = 1
            case 2
                threshold = 216;    %level = 2
            case 3
                threshold = 0;
        end
        for j = 1:compresslevel
            cDh = detcoef2('h',Cmod,S,j);
            [coeffrows coeffcolumns] = size(cDh);
            lower_cDh = find(abs(cDh(1:coeffrows,1:coeffcolumns)) <= threshold);
            cDh(lower_cDh) = 0;
            Cmod = modify_2D_wavelet_coeffs(Cmod,S,'h',j,cDh);
            cDv = detcoef2('v',Cmod,S,j);
            lower_cDv = find(abs(cDv(1:coeffrows,1:coeffcolumns)) <= threshold);
            cDv(lower_cDv) = 0;
            Cmod = modify_2D_wavelet_coeffs(Cmod,S,'v',j,cDv);
            cDd = detcoef2('d',Cmod,S,j);
            lower_cDd = find(abs(cDd(1:coeffrows,1:coeffcolumns)) <= threshold);
            cDd(lower_cDd) = 0;
            Cmod = modify_2D_wavelet_coeffs(Cmod,S,'d',j,cDd);
        end
           
    case {2,'LEVELWISE GLOBAL HARD THRESHOLDING OF H,V,D DETAIL COEFFS'}
        switch signaltype
            case 1                         % supply thresholds from level 1 to compresslevel
                thresholds = [37.22 46.27 0];
            case 2
                thresholds = [20 75 0];
            case 3
                thresholds = [0 0 0];
        end
        for j = 1:compresslevel
            cDh = detcoef2('h',Cmod,S,j);
            [coeffrows coeffcolumns] = size(cDh);
            lower_cDh = find(abs(cDh(1:coeffrows,1:coeffcolumns)) <= thresholds(j));
            cDh(lower_cDh) = 0;
            Cmod = modify_2D_wavelet_coeffs(Cmod,S,'h',j,cDh);
            cDv = detcoef2('v',Cmod,S,j);
            lower_cDv = find(abs(cDv(1:coeffrows,1:coeffcolumns)) <= thresholds(j));
            cDv(lower_cDv) = 0;
            Cmod = modify_2D_wavelet_coeffs(Cmod,S,'v',j,cDv);
            cDd = detcoef2('d',Cmod,S,j);
            lower_cDd = find(abs(cDd(1:coeffrows,1:coeffcolumns)) <= thresholds(j));
            cDd(lower_cDd) = 0;
            Cmod = modify_2D_wavelet_coeffs(Cmod,S,'d',j,cDd);
        end
end

Imod = waverec2(Cmod,S,wavelet);
%------------------------------------------------------------------------------------------------------------
% plot modified DFT/FWT coefficients at actual length
figure('Name',['ANALYSIS: compressed signal and decomposition in DWT/FWT coefficients for analysis with wavelet = '...
    wavelet , ' at level = ' num2str(compresslevel) '    (coefficients shown at correct dyadic size ~ border distortion due to FWT removed)'])
subplot(1,2,1)
image(wcodemat(Imod,nbcol));colormap(map);colorbar;axis image,title(['signal, size = [' num2str(rows) ' ' num2str(columns) ']'])
subplot(1,2,2)
j = 1:compresslevel;
lines_x_coord = [columns./2.^j+1 zeros(1,compresslevel); columns./2.^j+1 columns./2.^(j-1)];
lines_y_coord = [zeros(1,compresslevel) rows./2.^j+1; rows./2.^(j-1) rows./2.^j+1];
for j = compresslevel:-1:1
    if j == compresslevel
        cA = appcoef2(Cmod,S,wavelet,compresslevel);
        row_shrinkage = S(1,1) - rows/2^compresslevel;
        row_shrinkage_top = floor(row_shrinkage/2);
        column_shrinkage = S(1,2) - columns/2^compresslevel;
        column_shrinkage_left = floor(column_shrinkage/2);
        cA_central = cA((row_shrinkage_top + 1):(row_shrinkage_top + rows/2^compresslevel),(column_shrinkage_left + 1):(column_shrinkage_left + columns/2^compresslevel)); 
    end
    cDh = detcoef2('h',Cmod,S,j);
    row_shrinkage = S(2 + compresslevel - j,1) - rows/2^j;
    row_shrinkage_top = floor(row_shrinkage/2);
    column_shrinkage = S(2 + compresslevel - j,2) - columns/2^j;
    column_shrinkage_left = floor(column_shrinkage/2);
    cDh_central = cDh((row_shrinkage_top + 1):(row_shrinkage_top + rows/2^j),(column_shrinkage_left + 1):(column_shrinkage_left + columns/2^j)); 
    cDv = detcoef2('v',Cmod,S,j);
    cDv_central = cDv((row_shrinkage_top + 1):(row_shrinkage_top + rows/2^j),(column_shrinkage_left + 1):(column_shrinkage_left + columns/2^j)); 
    cDd = detcoef2('d',Cmod,S,j);
    cDd_central = cDd((row_shrinkage_top + 1):(row_shrinkage_top + rows/2^j),(column_shrinkage_left + 1):(column_shrinkage_left + columns/2^j)); 
    if j == compresslevel
        Icoeffs = [wcodemat(cA_central,nbcol) wcodemat(cDh_central,nbcol); wcodemat(cDv_central,nbcol) wcodemat(cDd_central,nbcol)];
    else
        Icoeffs = [Icoeffs wcodemat(cDh_central,nbcol); wcodemat(cDv_central,nbcol) wcodemat(cDd_central,nbcol)];
    end
end
image(wcodemat(Icoeffs,nbcol));colormap(map);colorbar;axis image,title(['cA: approx. coeffs [' num2str(rows) '/2^l^e^v^e^l ' num2str(columns) '/2^l^e^v^e^l]                  cDh: horisontal detail coeffs [' num2str(rows) '/2^l^e^v^e^l ' num2str(columns) '/2^l^e^v^e^l]']),...
    xlabel(['cDv: vertical detail coeffs [' num2str(rows) '/2^l^e^v^e^l ' num2str(columns) '/2^l^e^v^e^l]           cDd: diagonal detail coeffs [' num2str(rows) '/2^l^e^v^e^l ' num2str(columns) '/2^l^e^v^e^l]']),line(lines_x_coord,lines_y_coord,'LineWidth',2,'Color','b')
%------------------------------------------------------------------------------------------------------------
% plot modified DFT/FWT coefficients expanded to signal length
figure('Name',['ANALYSIS: compressed signal and decomposition in DWT/FWT coefficients for analysis with wavelet = '...
    wavelet , '    (complete coefficients including border distortion due to FWT, shown at signal size)'])
subplot(compresslevel+1,4,1)
image(wcodemat(Imod,nbcol));colormap(map);colorbar;axis image,title(['signal, size = [' num2str(rows) ' ' num2str(columns) ']'])
for j = 1:compresslevel
    subplot(compresslevel+1,4,j*4+1)
    cA = appcoef2(Cmod,S,wavelet,j);
    image(wcodemat(cA,nbcol));colormap(map);colorbar;axis image,ylabel(['level = ' num2str(j)']),title(['cA: approx. coeffs [' num2str(size(cA,1)) ' ' num2str(size(cA,2)) ']'])
    subplot(compresslevel+1,4,j*4+2)
    cDh = detcoef2('h',Cmod,S,j);
    image(wcodemat(cDh,nbcol));colormap(map);colorbar;axis image,title(['cDh: horisontal detail coeffs [' num2str(size(cDh,1)) ' ' num2str(size(cDh,2)) ']'])
    subplot(compresslevel+1,4,j*4+3)
    cDv = detcoef2('v',Cmod,S,j);
    image(wcodemat(cDv,nbcol));colormap(map);colorbar;axis image,title(['cDv: vertical detail coeffs [' num2str(size(cDv,1)) ' ' num2str(size(cDv,2)) ']'])
    subplot(compresslevel+1,4,j*4+4)
    cDd = detcoef2('d',Cmod,S,j);
    image(wcodemat(cDd,nbcol));colormap(map);colorbar;axis image,title(['cDd: diagonal detail coeffs [' num2str(size(cDd,1)) ' ' num2str(size(cDd,2)) ']'])
end
%-----------------------------------------------------------------------------------------------------------
% plot modified reconstructed approximations and details
figure('Name',['SYNTHESIS: compressed signal and contributes from individually approximations and details reconstructed with wavelet = '...
    wavelet ' and truncated to correct signal size ~ border distortion due to IFWT removed'])
subplot(compresslevel+1,4,1)
image(wcodemat(Imod,nbcol));colormap(map);colorbar;axis image,title(['signal, size = [' num2str(rows) ' ' num2str(columns) ']'])
for j = 1:compresslevel
    subplot(compresslevel+1,4,j*4+1)
    A = wrcoef2('a',Cmod,S,wavelet,j);
    image(wcodemat(A,nbcol));colormap(map);colorbar;axis image,ylabel(['level = ' num2str(j)']),title(['A: approximation [' num2str(rows) ' ' num2str(columns) ']'])
    subplot(compresslevel+1,4,j*4+2)
    Dh = wrcoef2('h',Cmod,S,wavelet,j);
    image(wcodemat(Dh,nbcol));colormap(map);colorbar;axis image,title(['Dh: horisontal detail [' num2str(rows) ' ' num2str(columns) ']'])
    subplot(compresslevel+1,4,j*4+3)
    Dv = wrcoef2('v',Cmod,S,wavelet,j);
    image(wcodemat(Dv,nbcol));colormap(map);colorbar;axis image,title(['Dv: vertical detail [' num2str(rows) ' ' num2str(columns) ']'])
    subplot(compresslevel+1,4,j*4+4)
    Dd = wrcoef2('d',Cmod,S,wavelet,j);
    image(wcodemat(Dd,nbcol));colormap(map);colorbar;axis image,title(['Dd: diagonal detail [' num2str(rows) ' ' num2str(columns) ']'])
end
% ---------------------------------------------------------------------------------------------
% calculate signal energy-ratio, coefficients-norm-recovery and compression-degree for current decomposition-compression-reconstruction scheme
cA = appcoef2(C,S,wavelet,compresslevel);
coeffs_norm = norm(cA)^2;
number_of_coeffs = size(cA,1)*size(cA,2);
for j = 1:compresslevel
    cDh = detcoef2('h',C,S,j);
    coeffs_norm = coeffs_norm + norm(cDh)^2;
    number_of_coeffs = number_of_coeffs + size(cDh,1)*size(cDh,2);
    cDv = detcoef2('v',C,S,j);
    coeffs_norm = coeffs_norm + norm(cDv)^2;
    number_of_coeffs = number_of_coeffs + size(cDv,1)*size(cDv,2);
    cDd = detcoef2('d',C,S,j);
    coeffs_norm = coeffs_norm + norm(cDd)^2;
    number_of_coeffs = number_of_coeffs + size(cDd,1)*size(cDd,2);
end
mod_coeffs_norm = norm(cA)^2;
number_of_zero_coeffs = 0;
for j = 1:compresslevel
    cDh = detcoef2('h',Cmod,S,j);
    mod_coeffs_norm = mod_coeffs_norm + norm(cDh)^2;
    number_of_zero_coeffs = number_of_zero_coeffs + length(find(cDh == 0));
    cDv = detcoef2('v',Cmod,S,j);
    mod_coeffs_norm = mod_coeffs_norm + norm(cDv)^2;
    number_of_zero_coeffs = number_of_zero_coeffs + length(find(cDv == 0));
    cDd = detcoef2('d',Cmod,S,j);
    mod_coeffs_norm = mod_coeffs_norm + norm(cDd)^2;
    number_of_zero_coeffs = number_of_zero_coeffs + length(find(cDd == 0));
end
energy_ratio = 100*(norm(Imod)/norm(I))^2;
coeffs_norm_recovery = 100*mod_coeffs_norm/coeffs_norm;
zero_coeffs = 100*number_of_zero_coeffs/number_of_coeffs;
compressfactor = number_of_coeffs/(number_of_coeffs - number_of_zero_coeffs); 
%--------------------------------------------------------------------------------------------------------------
%plot image and compressed image
figure('Name',['Original image and compressed image'])
subplot(1,2,1)
image(wcodemat(I,nbcol));colormap(map);colorbar;axis image,title(['original signal, size = [' num2str(rows) ' ' num2str(columns) ']'])
subplot(1,2,2)
image(wcodemat(Imod,nbcol));colormap(map);colorbar;axis image,title(['compressed signal, size = [' num2str(rows) ' ' num2str(columns) ']']), ...
        text(50,rows+50,['Energy ratio = ' num2str(energy_ratio,4) '%']),text(50,rows+75,['Coeffs norm recovery = ' num2str(coeffs_norm_recovery,4) '%']),...
        text(50,rows+100,['Zero coeffs = ' num2str(zero_coeffs,4) '%']),text(50,rows+125,['Compressfactor = ' num2str(compressfactor,4)])
%--------------------------------------------------------------------------------------------------------------

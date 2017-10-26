% illustration of 1D and 2D FWT boundary distortion when dealing with
% finite-length (practical) signals

% boundary distortion is analyzed as a function of
% - signal length (size)   (even or odd)
% - wavelet                (orthogonal or biorthogonal)
% - QMF filter length      (even or odd)
% - signal extension mode  (dwtmode)

clear;clc;close all;

% ----------------------------------------------------------------------------------------------
% 1D analysis
% ----------------------------------------------------------------------------------------------

% ---------------------------------------------------------------------------------------------
% define signal
N = 1024;
n = 0:N-1;
fs = 8000;
Ts = 1/fs;
t = (0:N-1)*Ts;
A1 = 1;
f1 = fs/400;
x = A1*sin(2*pi*f1*t+pi/6);
sigma = 0.1;
s = x + sigma*randn(1,N);
signallength = length(s);
%-----------------------------------------------------------------------------------------------
% No DSP on coefficients: Perfect reconstruction for all extension modes
% (however different coefficients)
wavelet = 'db6';
maxlevel = 5;
dwtmode('sym','nodisp');
[C,L] = wavedec(s, maxlevel, wavelet);
cD_sym = detcoef(C,L,3);
M = length(cD_sym);
m = 0:M-1;
s_sym = wrcoef('a',C,L,wavelet,0);
dwtmode('symw','nodisp');
[C,L] = wavedec(s, maxlevel, wavelet);
cD_symw = detcoef(C,L,3);
s_symw = wrcoef('a',C,L,wavelet,0);
dwtmode('zpd','nodisp');
[C,L] = wavedec(s, maxlevel, wavelet);
cD_zpd = detcoef(C,L,3);
s_zpd = wrcoef('a',C,L,wavelet,0);
dwtmode('ppd','nodisp');
[C,L] = wavedec(s, maxlevel, wavelet);
cD_ppd = detcoef(C,L,3);
s_per = wrcoef('a',C,L,wavelet,0);
figure('Name',['Perfect reconstruction (no DSP) for any of the FWT / IFWT extension modes  (but different boundary coefficients)        wavelet = ' wavelet ',  analysis-level = ' num2str(maxlevel)])
subplot(3,4,1)
plot(n,s,'k-'),ylabel('signal'),axis([0 signallength-1 min(s) max(s)]),grid,title('extension mode = sym 1/2')
subplot(3,4,5)
plot(n,s,'k-',n,s_sym,'r-'),ylabel('reconstructed signal'),axis([0 signallength-1 min(s) max(s)]),grid,legend('signal','reconstructed signal')
subplot(3,4,2)
plot(n,s,'k-'),axis([0 signallength-1 min(s) max(s)]),grid,title('extension mode = sym 1')
subplot(3,4,6)
plot(n,s,'k-',n,s_symw,'b-'),axis([0 signallength-1 min(s) max(s)]),grid,legend('signal','reconstructed signal')
subplot(3,4,3)
plot(n,s,'k-'),axis([0 signallength-1 min(s) max(s)]),grid,title('extension mode = zpd')
subplot(3,4,7)
plot(n,s,'k-',n,s_zpd,'g-'),axis([0 signallength-1 min(s) max(s)]),grid,legend('signal','reconstructed signal')
subplot(3,4,4)
plot(n,s,'k-'),axis([0 signallength-1 min(s) max(s)]),grid,title('extension mode = ppd')
subplot(3,4,8)
plot(n,s,'k-',n,s_per,'m-'),axis([0 signallength-1 min(s) max(s)]),grid,legend('signal','reconstructed signal')
subplot(3,4,[9 12])
plot(m,cD_sym,'r-','LineWidth',2),hold on,plot(m,cD_symw,'b-','LineWidth',2),plot(m,cD_zpd,'g-','LineWidth',2),plot(m,cD_ppd,'m-','LineWidth',2),hold off,ylabel('detail coeffs level 3'),axis([0 M-1 min(cD_sym) max(cD_sym)]),grid,...
    legend('sym 1/2','sym1','zpd','ppd','Location','North')
% ---------------------------------------------------------------------------------------------
% boundary distortion for orthogonal analysis (non-symmetric filters) and DSP on coefficients
wavelet = 'db10';
maxlevel = 7;
dwtmode('sym','nodisp');
[C,L] = wavedec(s, maxlevel, wavelet);
[s_sym,CXD,LXD] = wden(C,L,'minimaxi','s','sln',maxlevel,wavelet);
dwtmode('symw','nodisp');
[C,L] = wavedec(s, maxlevel, wavelet);
[s_symw,CXD,LXD] = wden(C,L,'minimaxi','s','sln',maxlevel,wavelet);
dwtmode('zpd','nodisp');
[C,L] = wavedec(s, maxlevel, wavelet);
[s_zpd,CXD,LXD] = wden(C,L,'minimaxi','s','sln',maxlevel,wavelet);
dwtmode('ppd','nodisp');
[C,L] = wavedec(s, maxlevel, wavelet);
[s_ppd,CXD,LXD] = wden(C,L,'minimaxi','s','sln',maxlevel,wavelet);
figure('Name',['Boundary distortion for ORTHOGONAL analysis (NON-symmetric filters) and DSP (noise reduction) on coefficients for different FWT / IFWT extension modes          wavelet = ' wavelet ',  analysis-level = ' num2str(maxlevel)])
subplot(3,4,1)
plot(n,s,'k-'),ylabel('noisy signal'),axis([0 signallength-1 min(s) max(s)]),grid,title('extension mode = sym 1/2')
subplot(3,4,5)
plot(n,x,'k-',n,s_sym,'r-'),ylabel('processed signal (noise reduction)'),axis([0 signallength-1 min(s) max(s)]),grid,legend('signal','processed signal','Location','South')
subplot(3,4,9)
plot(n,x,'k-',n,s_sym,'r-','LineWidth',2),ylabel('zoom of left boundary of processed signal'),axis([0 150 min(s) max(s)]),grid,legend('signal','processed signal','Location','South')
subplot(3,4,2)
plot(n,s,'k-'),axis([0 signallength-1 min(s) max(s)]),grid,title('extension mode = sym 1')
subplot(3,4,6)
plot(n,x,'k-',n,s_symw,'b-'),axis([0 signallength-1 min(s) max(s)]),grid,legend('signal','processed signal','Location','South')
subplot(3,4,10)
plot(n,x,'k-',n,s_symw,'b-','LineWidth',2),axis([0 150 min(s) max(s)]),grid,legend('signal','processed signal','Location','South')
subplot(3,4,3)
plot(n,s,'k-'),axis([0 signallength-1 min(s) max(s)]),grid,title('extension mode = zpd')
subplot(3,4,7)
plot(n,x,'k-',n,s_zpd,'g-'),axis([0 signallength-1 min(s) max(s)]),grid,legend('signal','processed signal','Location','South')
subplot(3,4,11)
plot(n,x,'k-',n,s_zpd,'g-','LineWidth',2),axis([0 150 min(s) max(s)]),grid,legend('signal','processed signal','Location','South')
subplot(3,4,4)
plot(n,s,'k-'),axis([0 signallength-1 min(s) max(s)]),grid,title('extension mode = ppd')
subplot(3,4,8)
plot(n,x,'k-',n,s_ppd,'m-'),axis([0 signallength-1 min(s) max(s)]),grid,legend('signal','processed signal','Location','South')
subplot(3,4,12)
plot(n,x,'k-',n,s_ppd,'m-','LineWidth',2),axis([0 150 min(s) max(s)]),grid,legend('signal','processed signal','Location','South')
% ---------------------------------------------------------------------------------------------
% boundary distortion for biorthogonal analysis (even-length symmetric filters) and DSP on coefficients
wavelet = 'bior3.9';
maxlevel = 5;
dwtmode('sym','nodisp');
[C,L] = wavedec(s, maxlevel, wavelet);
[s_sym,CXD,LXD] = wden(C,L,'minimaxi','s','sln',maxlevel,wavelet);
dwtmode('symw','nodisp');
[C,L] = wavedec(s, maxlevel, wavelet);
[s_symw,CXD,LXD] = wden(C,L,'minimaxi','s','sln',maxlevel,wavelet);
dwtmode('zpd','nodisp');
[C,L] = wavedec(s, maxlevel, wavelet);
[s_zpd,CXD,LXD] = wden(C,L,'minimaxi','s','sln',maxlevel,wavelet);
dwtmode('ppd','nodisp');
[C,L] = wavedec(s, maxlevel, wavelet);
[s_ppd,CXD,LXD] = wden(C,L,'minimaxi','s','sln',maxlevel,wavelet);
figure('Name',['Boundary distortion for BI-ORTHOGONAL analysis (SYMMETRIC / ANTI-SYMM. filters) and DSP (noise reduction) on coeffs for different FWT / IFWT extension modes       wavelet = ' wavelet ',  level = ' num2str(maxlevel)])
subplot(3,4,1)
plot(n,s,'k-'),ylabel('noisy signal'),axis([0 signallength-1 min(s) max(s)]),grid,title('extension mode = sym 1/2')
subplot(3,4,5)
plot(n,x,'k-',n,s_sym,'r-'),ylabel('processed signal (noise reduction)'),axis([0 signallength-1 min(s) max(s)]),grid,legend('signal','processed signal','Location','South')
subplot(3,4,9)
plot(n,x,'k-',n,s_sym,'r-','LineWidth',2),ylabel('zoom of left boundary of processed signal'),axis([0 150 min(s) max(s)]),grid,legend('signal','processed signal','Location','South'),...
    text(0,-1.8,'For biorthogonal analysis with EVEN-LENGTH FILTERS (bior3.9: 20,4,4,20) minimum boundary distortion is theoretically achieved for EXTENSION MODE "sym 1/2" (symmetric coefficients at all levels)'),...
    text(425,-2,'(in theory only for even-length signals - but of limited practical importance if signallength >> 1)')
subplot(3,4,2)
plot(n,s,'k-'),axis([0 signallength-1 min(s) max(s)]),grid,title('extension mode = sym 1')
subplot(3,4,6)
plot(n,x,'k-',n,s_symw,'b-'),axis([0 signallength-1 min(s) max(s)]),grid,legend('signal','processed signal','Location','South')
subplot(3,4,10)
plot(n,x,'k-',n,s_symw,'b-','LineWidth',2),axis([0 150 min(s) max(s)]),grid,legend('signal','processed signal','Location','South')
subplot(3,4,3)
plot(n,s,'k-'),axis([0 signallength-1 min(s) max(s)]),grid,title('extension mode = zpd')
subplot(3,4,7)
plot(n,x,'k-',n,s_zpd,'g-'),axis([0 signallength-1 min(s) max(s)]),grid,legend('signal','processed signal','Location','South')
subplot(3,4,11)
plot(n,x,'k-',n,s_zpd,'g-','LineWidth',2),axis([0 150 min(s) max(s)]),grid,legend('signal','processed signal','Location','South')
subplot(3,4,4)
plot(n,s,'k-'),axis([0 signallength-1 min(s) max(s)]),grid,title('extension mode = ppd')
subplot(3,4,8)
plot(n,x,'k-',n,s_ppd,'m-'),axis([0 signallength-1 min(s) max(s)]),grid,legend('signal','processed signal','Location','South')
subplot(3,4,12)
plot(n,x,'k-',n,s_ppd,'m-','LineWidth',2),axis([0 150 min(s) max(s)]),grid,legend('signal','processed signal','Location','South')
%--------------------------------------------------------------------------------------------------------------
% boundary distortion for biorthogonal analysis (odd-length symmetric filters) and DSP on coefficients
wavelet = 'bior6.8';
maxlevel = 5;
dwtmode('sym','nodisp');
[C,L] = wavedec(s, maxlevel, wavelet);
[s_sym,CXD,LXD] = wden(C,L,'minimaxi','s','sln',maxlevel,wavelet);
dwtmode('symw','nodisp');
[C,L] = wavedec(s, maxlevel, wavelet);
[s_symw,CXD,LXD] = wden(C,L,'minimaxi','s','sln',maxlevel,wavelet);
dwtmode('zpd','nodisp');
[C,L] = wavedec(s, maxlevel, wavelet);
[s_zpd,CXD,LXD] = wden(C,L,'minimaxi','s','sln',maxlevel,wavelet);
dwtmode('ppd','nodisp');
[C,L] = wavedec(s, maxlevel, wavelet);
[s_ppd,CXD,LXD] = wden(C,L,'minimaxi','s','sln',maxlevel,wavelet);
figure('Name',['Boundary distortion for BI-ORTHOGONAL analysis (SYMMETRIC / ANTI-SYMM. filters) and DSP (noise reduction) on coeffs for different FWT / IFWT extension modes       wavelet = ' wavelet ',  level = ' num2str(maxlevel)])
subplot(3,4,1)
plot(n,s,'k-'),ylabel('noisy signal'),axis([0 signallength-1 min(s) max(s)]),grid,title('extension mode = sym 1/2')
subplot(3,4,5)
plot(n,x,'k-',n,s_sym,'r-'),ylabel('processed signal (noise reduction)'),axis([0 signallength-1 min(s) max(s)]),grid,legend('signal','processed signal','Location','South')
subplot(3,4,9)
plot(n,x,'k-',n,s_sym,'r-','LineWidth',2),ylabel('zoom of left boundary of processed signal'),axis([0 150 min(s) max(s)]),grid,legend('signal','processed signal','Location','South'),...
    text(0,-1.8,'For biorthogonal analysis with ODD-LENGTH FILTERS (bior6.8: 17,11,11,17) minimum boundary distortion is theoretically achieved for EXTENSION MODE "sym 1" (symmetric coefficients at all levels)'),...
    text(425,-2,'(in theory only for odd-length signals - but of limited practical importance if signallength >> 1)')
subplot(3,4,2)
plot(n,s,'k-'),axis([0 signallength-1 min(s) max(s)]),grid,title('extension mode = sym 1')
subplot(3,4,6)
plot(n,x,'k-',n,s_symw,'b-'),axis([0 signallength-1 min(s) max(s)]),grid,legend('signal','processed signal','Location','South')
subplot(3,4,10)
plot(n,x,'k-',n,s_symw,'b-','LineWidth',2),axis([0 150 min(s) max(s)]),grid,legend('signal','processed signal','Location','South')
subplot(3,4,3)
plot(n,s,'k-'),axis([0 signallength-1 min(s) max(s)]),grid,title('extension mode = zpd')
subplot(3,4,7)
plot(n,x,'k-',n,s_zpd,'g-'),axis([0 signallength-1 min(s) max(s)]),grid,legend('signal','processed signal','Location','South')
subplot(3,4,11)
plot(n,x,'k-',n,s_zpd,'g-','LineWidth',2),axis([0 150 min(s) max(s)]),grid,legend('signal','processed signal','Location','South')
subplot(3,4,4)
plot(n,s,'k-'),axis([0 signallength-1 min(s) max(s)]),grid,title('extension mode = ppd')
subplot(3,4,8)
plot(n,x,'k-',n,s_ppd,'m-'),axis([0 signallength-1 min(s) max(s)]),grid,legend('signal','processed signal','Location','South')
subplot(3,4,12)
plot(n,x,'k-',n,s_ppd,'m-','LineWidth',2),axis([0 150 min(s) max(s)]),grid,legend('signal','processed signal','Location','South')
% ---------------------------------------------------------------------------------------------


% ----------------------------------------------------------------------------------------------
% 2D analysis
% ----------------------------------------------------------------------------------------------

%-----------------------------------------------------------------------------------------------------
% define or load signal                           
signaltype = 1;
switch signaltype
    case 1
        felt = ones(32,32);
        A = [100*felt 150*felt; 150*felt 100*felt];
        Iorig = repmat(A,4);
        [rows columns] = size(Iorig);
        nbcol = 256;
        colormap(gray(nbcol));
        map = colormap; 
        
    case 2
        load wbarb;  % --> X = wbarb; map; indexed image
        Iorig = X;
        [rows columns] = size(Iorig);
        nbcol = size(map,1);       
end
%---------------------------------------------------------------------------------------------------------
% add noise
sigma = 30;
noise = sigma*randn(rows,columns);
I = Iorig + noise;
%---------------------------------------------------------------------------------------------------------
% No DSP on coefficients: Perfect reconstruction for all extension modes
% (however different coefficients)
wavelet = 'db10';
maxlevel = 3;
dwtmode('sym','nodisp');
[C,S] = wavedec2(Iorig,maxlevel,wavelet);
cDh2_sym = detcoef2('h',C,S,2);
I_sym = wrcoef2('a',C,S,wavelet,0);
dwtmode('symw','nodisp');
[C,S] = wavedec2(Iorig,maxlevel,wavelet);
cDh2_symw = detcoef2('h',C,S,2);
I_symw = wrcoef2('a',C,S,wavelet,0);
dwtmode('zpd','nodisp');
[C,S] = wavedec2(Iorig,maxlevel,wavelet);
cDh2_zpd = detcoef2('h',C,S,2);
I_zpd = wrcoef2('a',C,S,wavelet,0);
dwtmode('ppd','nodisp');
[C,S] = wavedec2(Iorig,maxlevel,wavelet);
cDh2_ppd = detcoef2('h',C,S,2);
I_ppd = wrcoef2('a',C,S,wavelet,0);
figure('Name',['Perfect reconstruction (no DSP) for any of the FWT / IFWT extension modes  (but different boundary coefficients)        wavelet = ' wavelet ',  analysis-level = ' num2str(maxlevel)])
subplot(3,4,1)
image(wcodemat(Iorig,nbcol));colormap(map);colorbar;axis image,ylabel('signal'),title('extension mode = sym 1/2')
subplot(3,4,5)
image(wcodemat(I_sym,nbcol));colormap(map);colorbar;axis image,ylabel('reconstructed signal'),text(-100,-50,['max reconstruction error = ' num2str(max(max(abs(Iorig-I_sym))))])
subplot(3,4,9)
image(wcodemat(cDh2_sym,nbcol));colormap(map);colorbar;axis image,ylabel(['horisontal detail coeffs level 2'])
subplot(3,4,2)
image(wcodemat(Iorig,nbcol));colormap(map);colorbar;axis image,title('extension mode = sym 1')
subplot(3,4,6)
image(wcodemat(I_symw,nbcol));colormap(map);colorbar;axis image,text(-100,-50,['max reconstruction error = ' num2str(max(max(abs(Iorig-I_symw))))])
subplot(3,4,10)
image(wcodemat(cDh2_symw,nbcol));colormap(map);colorbar;axis image
subplot(3,4,3)
image(wcodemat(Iorig,nbcol));colormap(map);colorbar;axis image,title('extension mode = zpd')
subplot(3,4,7)
image(wcodemat(I_zpd,nbcol));colormap(map);colorbar;axis image,text(-100,-50,['max reconstruction error = ' num2str(max(max(abs(Iorig-I_zpd))))])
subplot(3,4,11)
image(wcodemat(cDh2_zpd,nbcol));colormap(map);colorbar;axis image
subplot(3,4,4)
image(wcodemat(Iorig,nbcol));colormap(map);colorbar;axis image,title('extension mode = ppd')
subplot(3,4,8)
image(wcodemat(I_ppd,nbcol));colormap(map);colorbar;axis image,text(-100,-50,['max reconstruction error = ' num2str(max(max(abs(Iorig-I_ppd))))])
subplot(3,4,12)
image(wcodemat(cDh2_ppd,nbcol));colormap(map);colorbar;axis image
% ---------------------------------------------------------------------------------------------
% boundary distortion for orthogonal analysis (non-symmetric filters) and
% DSP on coefficients
wavelet = 'db10';
maxlevel = 3;
dwtmode('sym','nodisp');
[thr,sorh,keepapp] = ddencmp('den','wv',I) ;
I_sym = wdencmp('gbl',I,wavelet,maxlevel,thr,sorh,keepapp);
zoom_rows = 1:64;
zoom_columns = 193:256;
I_sym_part = I_sym(zoom_rows,zoom_columns);
dwtmode('symw','nodisp');
I_symw = wdencmp('gbl',I,wavelet,maxlevel,thr,sorh,keepapp);
I_symw_part = I_symw(zoom_rows,zoom_columns);
dwtmode('zpd','nodisp');
I_zpd = wdencmp('gbl',I,wavelet,maxlevel,thr,sorh,keepapp);
I_zpd_part = I_zpd(zoom_rows,zoom_columns);
dwtmode('ppd','nodisp');
I_ppd = wdencmp('gbl',I,wavelet,maxlevel,thr,sorh,keepapp);
I_ppd_part = I_ppd(zoom_rows,zoom_columns);
figure('Name',['Boundary distortion for ORTHOGONAL analysis (NON-symmetric filters) and DSP (noise reduction) on coefficients for different FWT / IFWT extension modes          wavelet = ' wavelet ',  analysis-level = ' num2str(maxlevel)])
subplot(3,4,1)
image(wcodemat(I,nbcol));colormap(map);colorbar;axis image,ylabel('noisy signal'),title('extension mode = sym 1/2')
subplot(3,4,5)
image(wcodemat(I_sym,nbcol));colormap(map);colorbar;axis image,ylabel('processed signal (noise reduction)')
subplot(3,4,9)
image(wcodemat(I_sym_part,nbcol));colormap(map);colorbar;axis square,axis off,...
    title(['rows = ' num2str(min(zoom_rows)) ':' num2str(max(zoom_rows)) ', columns = ' num2str(min(zoom_columns)) ':' num2str(max(zoom_columns))]) 
subplot(3,4,2)
image(wcodemat(I,nbcol));colormap(map);colorbar;axis image,title('extension mode = sym 1')
subplot(3,4,6)
image(wcodemat(I_symw,nbcol));colormap(map);colorbar;axis image
subplot(3,4,10)
image(wcodemat(I_symw_part,nbcol));colormap(map);colorbar;axis square,axis off,...
    title(['rows = ' num2str(min(zoom_rows)) ':' num2str(max(zoom_rows)) ', columns = ' num2str(min(zoom_columns)) ':' num2str(max(zoom_columns))])
subplot(3,4,3)
image(wcodemat(I,nbcol));colormap(map);colorbar;axis image,title('extension mode = zpd')
subplot(3,4,7)
image(wcodemat(I_zpd,nbcol));colormap(map);colorbar;axis image
subplot(3,4,11)
image(wcodemat(I_zpd_part,nbcol));colormap(map);colorbar;axis square,axis off,...
    title(['rows = ' num2str(min(zoom_rows)) ':' num2str(max(zoom_rows)) ', columns = ' num2str(min(zoom_columns)) ':' num2str(max(zoom_columns))])
subplot(3,4,4)
image(wcodemat(I,nbcol));colormap(map);colorbar;axis image,title('extension mode = ppd')
subplot(3,4,8)
image(wcodemat(I_ppd,nbcol));colormap(map);colorbar;axis image
subplot(3,4,12)
image(wcodemat(I_ppd_part,nbcol));colormap(map);colorbar;axis square,axis off,...
    title(['rows = ' num2str(min(zoom_rows)) ':' num2str(max(zoom_rows)) ', columns = ' num2str(min(zoom_columns)) ':' num2str(max(zoom_columns))])
%--------------------------------------------------------------------------------------------------------------
% boundary distortion for biorthogonal analysis (even-length symmetric% filters) and DSP on coefficients
wavelet = 'bior3.9';
maxlevel = 3;
dwtmode('sym','nodisp');
[thr,sorh,keepapp] = ddencmp('den','wv',I) ;
I_sym = wdencmp('gbl',I,wavelet,maxlevel,thr,sorh,keepapp);
zoom_rows = 1:64;
zoom_columns = 193:256;
I_sym_part = I_sym(zoom_rows,zoom_columns);
dwtmode('symw','nodisp');
I_symw = wdencmp('gbl',I,wavelet,maxlevel,thr,sorh,keepapp);
I_symw_part = I_symw(zoom_rows,zoom_columns);
dwtmode('zpd','nodisp');
I_zpd = wdencmp('gbl',I,wavelet,maxlevel,thr,sorh,keepapp);
I_zpd_part = I_zpd(zoom_rows,zoom_columns);
dwtmode('ppd','nodisp');
I_ppd = wdencmp('gbl',I,wavelet,maxlevel,thr,sorh,keepapp);
I_ppd_part = I_ppd(zoom_rows,zoom_columns);
figure('Name',['Boundary distortion for BI-ORTHOGONAL analysis (SYMMETRIC / ANTI-SYMM. filters) and DSP (noise reduction) on coeffs for different FWT / IFWT extension modes       wavelet = ' wavelet ',  level = ' num2str(maxlevel)])
subplot(3,4,1)
image(wcodemat(I,nbcol));colormap(map);colorbar;axis image,ylabel('noisy signal'),title('extension mode = sym 1/2')
subplot(3,4,5)
image(wcodemat(I_sym,nbcol));colormap(map);colorbar;axis image,ylabel('processed signal (noise reduction)')
subplot(3,4,9)
image(wcodemat(I_sym_part,nbcol));colormap(map);colorbar;axis square,axis off,...
    title(['rows = ' num2str(min(zoom_rows)) ':' num2str(max(zoom_rows)) ', columns = ' num2str(min(zoom_columns)) ':' num2str(max(zoom_columns))]),...
    text(0,100,'For biorthogonal analysis with EVEN-LENGTH FILTERS (bior3.9: 20,4,4,20) minimum boundary distortion is theoretically achieved for EXTENSION MODE "sym 1/2" (symmetric coefficients at all levels)'),...
    text(350,110,'(in theory only for even-length signals - but of limited practical importance if signallength >> 1)')
subplot(3,4,2)
image(wcodemat(I,nbcol));colormap(map);colorbar;axis image,title('extension mode = sym 1')
subplot(3,4,6)
image(wcodemat(I_symw,nbcol));colormap(map);colorbar;axis image
subplot(3,4,10)
image(wcodemat(I_symw_part,nbcol));colormap(map);colorbar;axis square,axis off,...
    title(['rows = ' num2str(min(zoom_rows)) ':' num2str(max(zoom_rows)) ', columns = ' num2str(min(zoom_columns)) ':' num2str(max(zoom_columns))])
subplot(3,4,3)
image(wcodemat(I,nbcol));colormap(map);colorbar;axis image,title('extension mode = zpd')
subplot(3,4,7)
image(wcodemat(I_zpd,nbcol));colormap(map);colorbar;axis image
subplot(3,4,11)
image(wcodemat(I_zpd_part,nbcol));colormap(map);colorbar;axis square,axis off,...
    title(['rows = ' num2str(min(zoom_rows)) ':' num2str(max(zoom_rows)) ', columns = ' num2str(min(zoom_columns)) ':' num2str(max(zoom_columns))])
subplot(3,4,4)
image(wcodemat(I,nbcol));colormap(map);colorbar;axis image,title('extension mode = ppd')
subplot(3,4,8)
image(wcodemat(I_ppd,nbcol));colormap(map);colorbar;axis image
subplot(3,4,12)
image(wcodemat(I_ppd_part,nbcol));colormap(map);colorbar;axis square,axis off,...
    title(['rows = ' num2str(min(zoom_rows)) ':' num2str(max(zoom_rows)) ', columns = ' num2str(min(zoom_columns)) ':' num2str(max(zoom_columns))])
%--------------------------------------------------------------------------------------------------------------
% boundary distortion for biorthogonal analysis (odd-length symmetric% filters) and DSP on coefficients
wavelet = 'bior6.8';
maxlevel = 3;
dwtmode('sym','nodisp');
[thr,sorh,keepapp] = ddencmp('den','wv',I) ;
I_sym = wdencmp('gbl',I,wavelet,maxlevel,thr,sorh,keepapp);
zoom_rows = 1:64;
zoom_columns = 193:256;
I_sym_part = I_sym(zoom_rows,zoom_columns);
dwtmode('symw','nodisp');
I_symw = wdencmp('gbl',I,wavelet,maxlevel,thr,sorh,keepapp);
I_symw_part = I_symw(zoom_rows,zoom_columns);
dwtmode('zpd','nodisp');
I_zpd = wdencmp('gbl',I,wavelet,maxlevel,thr,sorh,keepapp);
I_zpd_part = I_zpd(zoom_rows,zoom_columns);
dwtmode('ppd','nodisp');
I_ppd = wdencmp('gbl',I,wavelet,maxlevel,thr,sorh,keepapp);
I_ppd_part = I_ppd(zoom_rows,zoom_columns);
figure('Name',['Boundary distortion for BI-ORTHOGONAL analysis (SYMMETRIC / ANTI-SYMM. filters) and DSP (noise reduction) on coeffs for different FWT / IFWT extension modes       wavelet = ' wavelet ',  level = ' num2str(maxlevel)])
subplot(3,4,1)
image(wcodemat(I,nbcol));colormap(map);colorbar;axis image,ylabel('noisy signal'),title('extension mode = sym 1/2')
subplot(3,4,5)
image(wcodemat(I_sym,nbcol));colormap(map);colorbar;axis image,ylabel('processed signal (noise reduction)')
subplot(3,4,9)
image(wcodemat(I_sym_part,nbcol));colormap(map);colorbar;axis square,axis off,...
    title(['rows = ' num2str(min(zoom_rows)) ':' num2str(max(zoom_rows)) ', columns = ' num2str(min(zoom_columns)) ':' num2str(max(zoom_columns))]),...
    text(0,100,'For biorthogonal analysis with ODD-LENGTH FILTERS (bior6.8: 17,11,11,17) minimum boundary distortion is theoretically achieved for EXTENSION MODE "sym 1" (symmetric coefficients at all levels)'),...
    text(350,110,'(in theory only for odd-length signals - but of limited practical importance if signallength >> 1)')
subplot(3,4,2)
image(wcodemat(I,nbcol));colormap(map);colorbar;axis image,title('extension mode = sym 1')
subplot(3,4,6)
image(wcodemat(I_symw,nbcol));colormap(map);colorbar;axis image
subplot(3,4,10)
image(wcodemat(I_symw_part,nbcol));colormap(map);colorbar;axis square,axis off,...
    title(['rows = ' num2str(min(zoom_rows)) ':' num2str(max(zoom_rows)) ', columns = ' num2str(min(zoom_columns)) ':' num2str(max(zoom_columns))])
subplot(3,4,3)
image(wcodemat(I,nbcol));colormap(map);colorbar;axis image,title('extension mode = zpd')
subplot(3,4,7)
image(wcodemat(I_zpd,nbcol));colormap(map);colorbar;axis image
subplot(3,4,11)
image(wcodemat(I_zpd_part,nbcol));colormap(map);colorbar;axis square,axis off,...
    title(['rows = ' num2str(min(zoom_rows)) ':' num2str(max(zoom_rows)) ', columns = ' num2str(min(zoom_columns)) ':' num2str(max(zoom_columns))])
subplot(3,4,4)
image(wcodemat(I,nbcol));colormap(map);colorbar;axis image,title('extension mode = ppd')
subplot(3,4,8)
image(wcodemat(I_ppd,nbcol));colormap(map);colorbar;axis image
subplot(3,4,12)
image(wcodemat(I_ppd_part,nbcol));colormap(map);colorbar;axis square,axis off,...
    title(['rows = ' num2str(min(zoom_rows)) ':' num2str(max(zoom_rows)) ', columns = ' num2str(min(zoom_columns)) ':' num2str(max(zoom_columns))])
%----------------------------------------------------------------------------------------------------------------
dwtmode('sym','nodisp');
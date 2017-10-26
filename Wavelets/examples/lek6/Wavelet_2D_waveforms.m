% illustration of 1D and 2D-separable wavelet functions

clear;clc;close all;

wavelet = 'db1';
iterations = 3;
[phi,psi,xval] = wavefun(wavelet,8);
[phiphi,phipsi,psiphi,psipsi,xyval] = wavefun2(wavelet,iterations);
[x y] = meshgrid(min(min(xyval)):(max(max(xyval))-min(min(xyval)))/(size(xyval,1)-1):max(max(xyval)));
figure('Name',['1D and 2D-separable tensor-product DWT/IDWT waveforms for orthogonal wavelet  "' wavelet '"'])
subplot(2,3,1)
plot(xval,phi,'b','LineWidth',2),grid,title('1D scaling function:   {\phi}(x)'),xlabel('x')
subplot(2,3,4)
plot(xval,psi,'r','LineWidth',2),grid,title('1D wavelet function:   {\psi}(x)'),xlabel('x')
subplot(2,3,2)
surf(x,y,phiphi),title('2D scaling function:     {\phi} (x,y) = {\phi}(x) {\phi}(y)')
subplot(2,3,3)
surf(x,y,phipsi),title('2D wavelet function 1:    {\psi}{\_}horisontal (x,y) = {\phi}(x) {\psi}(y)')
subplot(2,3,5)
surf(x,y,psiphi),title('2D wavelet function 2:    {\psi}{\_}vertical (x,y) = {\psi}(x) {\phi}(y)')
subplot(2,3,6)
surf(x,y,psipsi),title('2D wavelet function 3:    {\psi}{\_}diagonal (x,y) = {\psi}(x) {\psi}(y)')

wavelet = 'db2';
iterations = 3;
[phi,psi,xval] = wavefun(wavelet,8);
[phiphi,phipsi,psiphi,psipsi,xyval] = wavefun2(wavelet,iterations);
[x y] = meshgrid(min(min(xyval)):(max(max(xyval))-min(min(xyval)))/(size(xyval,1)-1):max(max(xyval)));
figure('Name',['1D and 2D-separable tensor-product DWT/IDWT waveforms for orthogonal wavelet  "' wavelet '"'])
subplot(2,3,1)
plot(xval,phi,'b','LineWidth',2),grid,title('1D scaling function:   {\phi}(x)'),xlabel('x')
subplot(2,3,4)
plot(xval,psi,'r','LineWidth',2),grid,title('1D wavelet function:   {\psi}(x)'),xlabel('x')
subplot(2,3,2)
surf(x,y,phiphi),title('2D scaling function:     {\phi} (x,y) = {\phi}(x) {\phi}(y)')
subplot(2,3,3)
surf(x,y,phipsi),title('2D wavelet function 1:    {\psi}{\_}horisontal (x,y) = {\phi}(x) {\psi}(y)')
subplot(2,3,5)
surf(x,y,psiphi),title('2D wavelet function 2:    {\psi}{\_}vertical (x,y) = {\psi}(x) {\phi}(y)')
subplot(2,3,6)
surf(x,y,psipsi),title('2D wavelet function 3:    {\psi}{\_}diagonal (x,y) = {\psi}(x) {\psi}(y)')

wavelet = 'db4';
iterations = 2;
[phi,psi,xval] = wavefun(wavelet,8);
[phiphi,phipsi,psiphi,psipsi,xyval] = wavefun2(wavelet,iterations);
[x y] = meshgrid(min(min(xyval)):(max(max(xyval))-min(min(xyval)))/(size(xyval,1)-1):max(max(xyval)));
figure('Name',['1D and 2D-separable tensor-product DWT/IDWT waveforms for orthogonal wavelet  "' wavelet '"'])
subplot(2,3,1)
plot(xval,phi,'b','LineWidth',2),grid,title('1D scaling function:   {\phi}(x)'),xlabel('x')
subplot(2,3,4)
plot(xval,psi,'r','LineWidth',2),grid,title('1D wavelet function:   {\psi}(x)'),xlabel('x')
subplot(2,3,2)
surf(x,y,phiphi),title('2D scaling function:     {\phi} (x,y) = {\phi}(x) {\phi}(y)')
subplot(2,3,3)
surf(x,y,phipsi),title('2D wavelet function 1:    {\psi}{\_}horisontal (x,y) = {\phi}(x) {\psi}(y)')
subplot(2,3,5)
surf(x,y,psiphi),title('2D wavelet function 2:    {\psi}{\_}vertical (x,y) = {\psi}(x) {\phi}(y)')
subplot(2,3,6)
surf(x,y,psipsi),title('2D wavelet function 3:    {\psi}{\_}diagonal (x,y) = {\psi}(x) {\psi}(y)')

wavelet = 'sym4';
iterations = 3;
[phi,psi,xval] = wavefun(wavelet,8);
[phiphi,phipsi,psiphi,psipsi,xyval] = wavefun2(wavelet,iterations);
[x y] = meshgrid(min(min(xyval)):(max(max(xyval))-min(min(xyval)))/(size(xyval,1)-1):max(max(xyval)));
figure('Name',['1D and 2D-separable tensor-product DWT/IDWT waveforms for orthogonal wavelet  "' wavelet '"'])
subplot(2,3,1)
plot(xval,phi,'b','LineWidth',2),grid,title('1D scaling function:   {\phi}(x)'),xlabel('x')
subplot(2,3,4)
plot(xval,psi,'r','LineWidth',2),grid,title('1D wavelet function:   {\psi}(x)'),xlabel('x')
subplot(2,3,2)
surf(x,y,phiphi),title('2D scaling function:     {\phi} (x,y) = {\phi}(x) {\phi}(y)')
subplot(2,3,3)
surf(x,y,phipsi),title('2D wavelet function 1:    {\psi}{\_}horisontal (x,y) = {\phi}(x) {\psi}(y)')
subplot(2,3,5)
surf(x,y,psiphi),title('2D wavelet function 2:    {\psi}{\_}vertical (x,y) = {\psi}(x) {\phi}(y)')
subplot(2,3,6)
surf(x,y,psipsi),title('2D wavelet function 3:    {\psi}{\_}diagonal (x,y) = {\psi}(x) {\psi}(y)')

wavelet = 'coif1';
iterations = 3;
[phi,psi,xval] = wavefun(wavelet,8);
[phiphi,phipsi,psiphi,psipsi,xyval] = wavefun2(wavelet,iterations);
[x y] = meshgrid(min(min(xyval)):(max(max(xyval))-min(min(xyval)))/(size(xyval,1)-1):max(max(xyval)));
figure('Name',['1D and 2D-separable tensor-product DWT/IDWT waveforms for orthogonal wavelet  "' wavelet '"'])
subplot(2,3,1)
plot(xval,phi,'b','LineWidth',2),grid,title('1D scaling function:   {\phi}(x)'),xlabel('x')
subplot(2,3,4)
plot(xval,psi,'r','LineWidth',2),grid,title('1D wavelet function:   {\psi}(x)'),xlabel('x')
subplot(2,3,2)
surf(x,y,phiphi),title('2D scaling function:     {\phi} (x,y) = {\phi}(x) {\phi}(y)')
subplot(2,3,3)
surf(x,y,phipsi),title('2D wavelet function 1:    {\psi}{\_}horisontal (x,y) = {\phi}(x) {\psi}(y)')
subplot(2,3,5)
surf(x,y,psiphi),title('2D wavelet function 2:    {\psi}{\_}vertical (x,y) = {\psi}(x) {\phi}(y)')
subplot(2,3,6)
surf(x,y,psipsi),title('2D wavelet function 3:    {\psi}{\_}diagonal (x,y) = {\psi}(x) {\psi}(y)')
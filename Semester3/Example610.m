clear;
%Example 6.10
%Frequensy and magnitude plot

% y(n) = 0.5*x(n)+0.5*x(n-1)

% Determine the frequency response

% y(z) = 0.5*x(z)+0.5x(z)z^-1
% H(z) = 0.5 + 0.5*z^(-1)

% H( e^(jw) ) = 0.5 + 0.5*e^(jwT)           //Substitue z = e^(jw)
% H( e^(jw) ) = 0.5 + cos(jw) + j*sin(jw)   //Applyed the qeuler formula

%|H( e^(jwT) )| =  sqrt(H( e^(jwT) ) * H( e^(jwT) )(conjugate))


figure (1)
[h,w] = freqz([1 -0.5], [1], 1024)
phi = 180*unwrap (angle(h) ) /pi;
subplot(2,1,1), plot(w,abs(h)),grid;xlabel('Frequensy radians');ylabel('Magnitude')
subplot(2,1,2), plot(w,phi),grid;xlabel('Frequensy radians');ylabel('Phase degrees')

figure (2)
[h,w] = freqz([1], [1 -0.5], 1024)
phi = 180*unwrap (angle(h) ) /pi;
subplot(2,1,1), plot(w,abs(h)),grid;xlabel('Frequensy radians');ylabel('Magnitude')
subplot(2,1,2), plot(w,phi),grid;xlabel('Frequensy radians');ylabel('Phase degrees')

figure (3)
[h,w] = freqz([0.5 0 -0.3], [1 -0.5 0.25], 1024)
phi = 180*unwrap (angle(h) ) /pi;
subplot(2,1,1), plot(w,abs(h)),grid;xlabel('Frequensy radians');ylabel('Magnitude')
subplot(2,1,2), plot(w,phi),grid;xlabel('Frequensy radians');ylabel('Phase degrees')

figure (4)
[h,w] = freqz([1 -0.9 0.81], [1 -0.6 0.36], 1024)
phi = 180*unwrap (angle(h) ) /pi;
subplot(2,1,1), plot(w,abs(h)),grid;xlabel('Frequensy radians');ylabel('Magnitude')
subplot(2,1,2), plot(w,phi),grid;xlabel('Frequensy radians');ylabel('Phase degrees')

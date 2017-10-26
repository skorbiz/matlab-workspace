clear;
%Example 6.1
%Compute Recursive function

%y(n) = 0.5y(n-2)+x(n-1)

%Initial conditions:
% y(-2) = 1
% y(-1) = 0
% x(-1) = -1

%Input
% x(n) = (0.5)^n*u(n)

y = zeros(1,10);       % setup vectpr to store y(n)
y = [1 0 y];           % Set initial conditions of y(-2) and y(-1)
n = 0:1:9;             % Compute time index
x = (0.5).^n;          % Compute 20 input samples of x(n)
x = [0 -1 x];          % set initial conditions of x(-2) = 0 and x(-1) = 1

for n =1:10
y(n +2) = 0.5*y(n)+x(n+1)
end

n = 0:1:9;
subplot(3,1,1);stem(n,x(3:12));ylabel('Input x(n)');xlabel('Sample number')
subplot(3,1,2);stem(n,y(3:12));ylabel('Output x(n)');xlabel('Sample number')
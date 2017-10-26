% Integrating Ordinary Differential equations
    clear all
global  x  xd  xdd

% Initial conditions
    x  = [1; -1];
    xd = [0;  1];
% Initialize the integration array with corrected I.C.'s
    u0 = [x; xd];
% Input time data
     dt = input('Enter time-step in seconds = ? ');
     t0 = 0;
     te = input('Enter final time = ? ');
     Tspan = [t0:dt:te];
% Begin integration
     [T, uT] = ode45('diffeq', Tspan, u0);
% The output T is a column array
% T(i) corresponds to a time mark
% uT(i,:) is a time record of uT' at T(i)

% Plot the results
    plot_sol (T,uT)

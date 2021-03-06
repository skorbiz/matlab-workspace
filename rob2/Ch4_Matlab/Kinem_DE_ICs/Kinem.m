% A Four-bar Mechanism
% Kinematic analysis by integration
% Correct initial conditions are computed
    clear all
    
global  l1 l2 l3
global  r1 r2 r3 r4
global  theta1 theta2 theta3
global  theta1d theta2d theta3d
global  theta1dd theta2dd theta3dd

% ----------------------------------------------------------------------
%                           Enter Data block
% ----------------------------------------------------------------------

% Constant coordinates, lengths, vectors
    l1 = 1;
    l2 = 3;
    l3 = 2.2;
    a  = 2;
    b  = 0.5;
% Estimates for variable coordinates
    theta1 = 90; % degrees
    theta2 = 30; % degrees
    theta3 = 75; % degrees

% ----------------------------------------------------------------------
%                 Initial Conditions Correction block
% ----------------------------------------------------------------------

    r4 = [a; b];
    theta1 = theta1*pi/180;
    theta2 = theta2*pi/180;
    theta3 = theta3*pi/180;
% Initial Condition correction for coordinate
    newton
% ' Corrected coordinates'
%    theta1
%    theta2
%    theta3
% Initial Condition correction for velocities
    vels
% ' Corrected velocities'    
%    theta1d
%    theta2d
%    theta3d

% ----------------------------------------------------------------------
%                           Integration block
% ----------------------------------------------------------------------

% Initialize the integration array with corrected I.C.'s
    u0 = [theta1 theta2 theta3 theta1d theta2d theta3d]';
    
% Input time data
     dt = input('Enter time-step in seconds = ? ');
     t_initial = 0;
     t_final = input('Enter final time = ? ');
     Tspan = [t_initial:dt:t_final];
     
% Begin integration
     [T, uT] = ode45 ('diffeq',Tspan,u0);

% The output T is a column vector; T(i) corresponds to a time mark
% uT(i,:) is a time record of uT' at T(i)

% ----------------------------------------------------------------------
%                       Post Processing block
% ----------------------------------------------------------------------

% Note: Accelerations at the integration time steps have been lost!!!
% We can recover (recompute) them if necessary

% We may report (values, plots, graphices, etc.) the contents of T and uT
  

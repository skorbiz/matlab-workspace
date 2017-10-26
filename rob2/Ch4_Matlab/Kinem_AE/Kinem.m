% A Four-bar Mechanism
% Kinematic analysis by solving algebraic equations
    clear all

global  l1 l2 l3
global  r1 r2 r3 r4
global  theta1 theta2 theta3
global  theta1d theta2d theta3d
global  theta1dd theta2dd theta3dd
global  theta1_0 w1_0 w1d_0
global  t_initial t_final dt t

% ----------------------------------------------------------------------
%                           Enter Data block
% ----------------------------------------------------------------------

% Constant coordinates, lengths, vectors
    l1 = 1;
    l2 = 3;
    l3 = 2.2;
    a  = 2;
    b  = 0.5;
% Driver data
    theta1_0 = pi/2; % rad
    w1_0    = 2*pi; % rad/sec
    w1d_0   = 0; % rad/sec^2
% Estimates for variable coordinates
    theta1 = 90; % degrees
    theta2 = 30; % degrees
    theta3 = 75; % degrees
% Input time data
    dt = input('Enter time-step in seconds = ? ');
    t_initial = 0;
    t_final = input('Enter final time = ? ');
% ----------------------------------------------------------------------

    r4 = [a; b];
    theta1 = theta1*pi/180;
    theta2 = theta2*pi/180;
    theta3 = theta3*pi/180;
% Determine number of time steps plus one
    n = (t_final - t_initial)/dt + 1;
% Initialize arrays for saving results
    T   = zeros (n,1);
    u   = zeros (n,3); ud  = zeros (n,3); udd = zeros (n,3);
    
% ----------------------------------------------------------------------
%                          Kinematic Analysis
% ----------------------------------------------------------------------

for i=1:n
    t = (i-1)*dt;
    T(i) = t;        
% Coordinate computation
    newton
    % save results
    u(i,1) = theta1; u(i,2) = theta2; u(i,3) = theta3;
% Velocity computation
    vels
    % save results
    ud(i,1) = theta1d; ud(i,2) = theta2d; ud(i,3) = theta3d;
% Acceleration analysis
    eqsmotion(t);
    % save results
    udd(i,1) = theta1dd; udd(i,2) = theta2dd; udd(i,3) = theta3dd;
end

% ----------------------------------------------------------------------

% We may report (values, plots, graphices, etc.) the contents of 
% T, u, ud and udd

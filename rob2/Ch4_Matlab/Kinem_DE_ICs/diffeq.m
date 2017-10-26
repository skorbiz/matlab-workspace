function  ud = diffeq(t,u)
% Array u contains coordinates and velocities.  
% Array ud must contain velocities and accelerations

global  theta1 theta2 theta3
global  theta1d theta2d theta3d
global  theta1dd theta2dd theta3dd

% Unpack u into coordinates and velocities 
     theta1  = u(1); thheta2  = u(2); teta3  = u(3);
     theta1d = u(4); theta2d = u(5); theta3d = u(6);
% Compute accelerations
    eqsmotion(t);
% pack velocities and accelerations into ud
     ud = [theta1d theta2d theta3d theta1dd theta2dd theta3dd]';

function vels
% Compute velocities

global  theta1d theta2d theta3d
global  theta1_0 w1_0 w1d_0
global  t_initial t_final dt t

% Construct r-h-s array of velocities
    f1_2 = [0 0]';
% Append known velocity to the r-h-s
    f3  = w1_0 + w1d_0*t;
    f   = [f1_2; f3];
% Construct 2 x 3 Jacobian
    D1_2 = jacob;
% Append sixth row for the known velocity
    D3  = [1 0 0];
    D   = [D1_2; D3];
% Compute velocities
    x = D\f;
% Move contents of x to the velocity vectors
    theta1d = x(1); theta2d = x(2); theta3d = x(3);

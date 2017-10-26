function rhs_A = rhsA
% Construct r-h-s array of acceleration constraints

global  r1 r2 r3 r4
global  theta1d theta2d theta3d

% Construct r-h-s array for acceleration constraints
     rhs_A = r1*theta1d^2 + r2*theta2d^2 - r3*theta3d^2;

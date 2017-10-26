function r_vectors
% Computes r vectors

global  l1 l2 l3
global  r1 r2 r3 r4
global  theta1 theta2 theta3

% Compute d and other vectors
    r1 = l1*[cos(theta1); sin(theta1)];
    r2 = l2*[cos(theta2); sin(theta2)];
    r3 = l3*[cos(theta3); sin(theta3)];
    
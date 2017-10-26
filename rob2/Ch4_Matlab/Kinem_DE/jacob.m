function D = jacob
% Construct Jacobian matrix

global  r1 r2 r3 r4

% Rotate vectors 90 degrees positively
     rot1 = rot(r1); rot2 = rot(r2); rot3 = rot(r3);
% Construct Jacobian
     D = [rot1 rot2 -rot3];

function f = cons
% Compute kinematic constraints

global  r1 r2 r3 r4

% Evaluate position constraints
    f = r1 + r2 - r3 - r4;



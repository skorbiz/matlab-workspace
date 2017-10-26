function xdd = eqs(x, xd, t)
% This routine computes accelerations
    
% Compute accelerations
    D = [1 0; 0 2];
    rhs = -[(2*xd(1) + 4*(x(1)-x(2))); ...
            (3*xd(2) + 4*(x(2)-x(1)))];
    xdd = D\rhs;
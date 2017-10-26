function  ud = diffeq(t,u)
% This routine receives array u from ODE45 and returns array ud (u_dot)
% to ODE45.  Array u contains coordinates and velocities.  Array ud must
% contain velocities and accelerations

% Unpack u into x and xd
     x  = u(1:2);
     xd = u(3:4);
% Compute xdd
     xdd = eqs(x, xd, t);
% pack r_d (velocities) and r_dd (accelerations) into ud
     ud = [xd; xdd];

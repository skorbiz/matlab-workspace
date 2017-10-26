% Solving nonlinear algebraic equations
% Method:   Newton-Raphson

global   l1 l2 l3 a b
global   theta1

% Constants
    l1 = 1.0; l2 = 3.0; l3 = 2.2; a = 2.2; b = 0.5;
% Known coordinate(s)
    theta1 = 90*pi/180;
% Initial estimates
    theta2 = 45*pi/180;
    theta3 = 80*pi/180;

% Move unknown coordinates to array x
    x = [theta2; theta3];

for n = 1:20
    % Compute 2 x 1 array of functions
        f = constraints (x);
    % condition for termination:
          normf = norm(f);
          if ( normf <= 1e-7 ) break; end;  
    % Construct 2 x 2 Jacobian
        D = jacobian(x);
    % Compute corrections
        delta_x = D\f;
    % Correct the coordinates
        x = x - delta_x;
end

% Report solution in degrees
    theta2 = x(1)*180/pi
    theta3 = x(2)*180/pi

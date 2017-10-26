function newton
% Newton-Raphson routine for solving nonlinear equations

global  theta1 theta2 theta3

% Move coordinates to array x
    x = [theta1; theta2; theta3];

for n = 1:20
    % Compute r vectors
        r_vectors;
    % Compute 2 x 1 array of functions
        f1_2 = cons;
    % Append known coordinate as a function
        f3  = theta1 - pi/2;
        f   = [f1_2; f3];
    % condition for termination:
          normf = norm(f);
          if ( normf <= 1e-7 ) break; end;  
    % Construct 2 x 3 Jacobian
        D1_2 = jacob;
    % Append sixth row for the known coordinate
        D3  = [1 0 0];
        D   = [D1_2; D3];
    % Compute corrections
        delta_x = D\f;
    % Correct coordinates
        x = x - delta_x;
    % Move contents of x to the coordinate vectors
        theta1 = x(1); theta2 = x(2); theta3 = x(3);
end

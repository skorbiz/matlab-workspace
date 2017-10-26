function  eqsmotion (t)
% Compute accelerations

global  theta1dd theta2dd theta3dd

% Update r vectors
    r_vectors
% Construct Jacobian
    D1_2 = jacob;
    D3   = [1 0 0];
    D    = [D1_2; D3];
% Construct r-h-s array of accelerations
    rhs1_2 = rhsA;
    rhs3   = 0;
    rhs    = [rhs1_2; rhs3];
% Compute accelerations
    x = D\rhs;

    theta1dd = x(1);
    theta2dd = x(2);
    theta3dd = x(3);
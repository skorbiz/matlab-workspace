% Solving a set of linear algebraic equations
% Method:   Matrix inversion

    A   = [ 3  1  -1
           -1  2   1
            2 -3   1 ];
    rhs = [ 2  6 -1 ]';

    x   = inv(A) * rhs
    
% Method:   Backslash command

    x   = A \ rhs
    
% Method:   L-U factorization
      
    [L, U, P] = lu(A)
        y = L \ (P * rhs)
        x = U \ y


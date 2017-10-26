function D = jacobian(x)
% This routine constructs the Jacobian matrix

global   l1 l2 l3 a b
global   theta1

theta2 = x(1);
theta3 = x(2);

d11 = -l2*sin(theta2); 
d12 =  l3*sin(theta3);
d21 =  l2*cos(theta2); 
d22 = -l3*cos(theta3);

D = [d11 d12; d21 d22];


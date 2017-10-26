function f = constraints(x)
% This routine evaluates the constraint equations

global   l1 l2 l3 a b
global   theta1

theta2 = x(1);
theta3 = x(2);

f1 = l1*cos(theta1) + l2*cos(theta2) - l3*cos(theta3) - a;
f2 = l1*sin(theta1) + l2*sin(theta2) - l3*sin(theta3) - b;

f = [f1 f2]';

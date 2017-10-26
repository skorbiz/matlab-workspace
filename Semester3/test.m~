%clear;
%syms z;
%f = (-0.25-1.25*z)*z/(0.25-1*z^2);
%s = iztrans(f);
%simplify(s)
%simple(s)
%f = (1.25*z*(z+0.2))/(z^2-0.25);
%[R,P,K] = residue([1.25 0.2*1.25 0],[1 0 -0.25])

%Input x(n) = direktDelta(n)
%Output y(n) = 0.25y(n-2) x(n)
%begyndelses betingelser y(-1) = y(-2) = 1
%ztrans(0.25*sym('y(n-2)')+sym('x(n)'))
%s = solve('yn=1/4+yn/(4*z^2)+1/(4*z)+1','yn')
%iztrans(s)

t = 0:pi/20:2*pi;
[x,y] = meshgrid(t);
z = (sin(x).^2)-(cos(y).^2);

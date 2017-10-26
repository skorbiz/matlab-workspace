
%% Problem 3.1
%

clc;clear;close all;
format long

max_iterationer=100;
epsilon = 1.0e-8;

x = 0;

disp('x in each ittteration:') 

for i = 1:max_iterationer
    
    f  = x^3+2*x-7;
    df = 3*x^2+2;
    xnew = x-f/df;
    if(abs(x-xnew) < epsilon)
        break;
    end
    x = xnew;
    disp(x)
    
end

disp('x equal to:')
disp(x)


%Plot the function
x = -10:0.01:10;

y = x.^3+2*x-7;
plot(x,y)
grid on % Turn on grid lines for this plot

%% Problem 3.2
%

clc;clear;close all;
format short

max_iterationer=100;
epsilon = 1.0e-8;

theta1 = pi*30/180;
q = [1; 0];

disp('d and theta2  in each ittteration:') 

for i = 1:max_iterationer
    phi1 = 10*cos(theta1) + 2*cos(q(2)) - q(1);
    phi2 = 10*sin(theta1) + 2*sin(q(2)) - 4;
    phi = [phi1; phi2];

    jac1 = [-1; 0];
    jac2 = [-2*sin(q(2)); 2*cos(q(2))];
    jaco  = [jac1 jac2];

    delta_q = -inv(jaco)*phi;

    q = q+delta_q;
    
    disp(q')

    if( norm(delta_q) < epsilon)
        break;
    end
end

disp('d and theta2 equal to:')
disp(q')

%% Problem 3.3
%

clc;clear;close all;
format short
%%
% Problem 1
% 

clc;clear;close all
format compact

d2r = pi/180;
 
vec_tilde = @(x) [0 -x(3) x(2); x(3) 0 -x(1); -x(2) x(1) 0]

%Transformation matrixes
A1 = @(x) [ 1 0 0; 0 cos(x*d2r) sin(x*d2r); 0 -sin(x*d2r) cos(x*d2r)];
A2 = @(x) [ cos(x*d2r) -sin(x*d2r) 0; sin(x*d2r)  cos(x*d2r) 0; 0 0 1];          %# An anonymous function

%The basis points in local coordinates
r1 = [0 0 0]'
S1q_local = [0 330 0]'
S2p_local = [300 0 0]'

%The global coordinates of the basis points (incl. point P):
S1q_global = r1 + A1(0)*S1q_local
S2p_global = S1q_global + A2(20)*S2p_local

%The angular velocities and accelerations of the actuators in local coordinates
w1_local = [4 0 0]'
w21_local = [0 0 1]'

w1_vel_local = [0 0 0]'
w21_vel_local = [0 0 0]'

%The angular velocities of the bodies
w1_global = A1(0)*w1_local
w2_global = w1_global + A1(0)*w21_local

%The angular acceleration of the bodies
w1_vel_global = [0 0 0]'
w2_vel_global = w1_vel_global + vec_tilde(w1_global)*A1(0)*w21_local + A1(0)*w21_vel_local

%The velocities of the basis points (incl P)
S1q_vel_local = [0 0 0]'
S2p_vel_local = [500 0 0]'
r1_vel = [0 0 0]';
r2_vel = r1_vel + vec_tilde(w1_global)*A1(0)*S1q_local + A1(0)*S1q_vel_local

S2p_vel_global = r2_vel + vec_tilde(w2_global)*A2(20)*S2p_local+A2(20)*S2p_vel_local

%The accelerations of the basis points (incl P)
S1q_acc_local = [0 0 0]'
r1_acc = [0 0 0]'
r2_acc = r1_acc + vec_tilde(w1_vel_global)*A1(0)*S1q_local ...
         + vec_tilde(w1_global)*vec_tilde(w1_global)*A1(0)*S1q_local ...
         + 2*vec_tilde(w1_global)*A1(0)*S1q_vel_local ...
         + A1(0)*S1q_acc_local

S2p_acc_global = r2_acc + vec_tilde(w2_vel_global)*A2(20)*S2p_local ...
         + vec_tilde(w2_global)*vec_tilde(w2_global)*A2(20)*S2p_local ...
         + 2*vec_tilde(w2_global)*A2(20)*S2p_vel_local ...
         + A2(20)*[0 0 0]'
     
     
%     



%%
% Problem 1
% 

clc;clear;close all
format compact

d2r = pi/180;
 
vec_tilde = @(x) [0 -x(3) x(2); x(3) 0 -x(1); -x(2) x(1) 0];

R_around_z = @(z) [ cos(z) -sin(z) 0; sin(z) cos(z) 0; 0 0 1]; 
R_around_y = @(y) [ cos(y) 0 sin(y); 0 1 0; -sin(y) 0 cos(y)]; 
R_around_x = @(x) [ 1 0 0; 0 cos(x) -sin(x); 0 sin(x) cos(x)]; 


%Transformation matrixes
A1_func = @(x) R_around_z(x*d2r);          %# An anonymous function
A2_func = @(x) R_around_z(x*d2r);

A1 = A1_func(60);
A2 = A2_func(30)';


%The basis points in local coordinates
S1q_local = [250 -100 0]'
S2p_local = [85 40 0]'

r1 = [0 0 0]'
r2 = r1 + A1*S1q_local

%The global coordinates of the basis points (incl. point P):
S1q_global = r1 + A1*S1q_local
S2p_global = r2 + A2*S2p_local

%The angular velocities and accelerations of the actuators in local coordinates
w1_local = [0 0 3]'
w21_local = [0 0 1]'

w1_vel_local = [0 0 0]'
w21_vel_local = [0 0 0]'

%The angular velocities of the bodies
w1_global = A1(0)*w1_local
w2_global = w1_global + A1(0)*w21_local

%The angular acceleration of the bodies
w1_vel_global = [0 0 0]'
w2_vel_global = w1_vel_global + vec_tilde(w1_global)*A1(0)*w21_local + A1(0)*w21_vel_local

%The velocities of the basis points (incl P)
S1q_vel_local = [0 0 0]'
S2p_vel_local = [500 0 0]'
r1_vel = [0 0 0]';
r2_vel = r1_vel + vec_tilde(w1_global)*A1(0)*S1q_local + A1(0)*S1q_vel_local

S2p_vel_global = r2_vel + vec_tilde(w2_global)*A2(20)*S2p_local+A2(20)*S2p_vel_local

%The accelerations of the basis points (incl P)
S1q_acc_local = [0 0 0]'
r1_acc = [0 0 0]'
r2_acc = r1_acc + vec_tilde(w1_vel_global)*A1(0)*S1q_local ...
         + vec_tilde(w1_global)*vec_tilde(w1_global)*A1(0)*S1q_local ...
         + 2*vec_tilde(w1_global)*A1(0)*S1q_vel_local ...
         + A1(0)*S1q_acc_local

S2p_acc_global = r2_acc + vec_tilde(w2_vel_global)*A2(20)*S2p_local ...
         + vec_tilde(w2_global)*vec_tilde(w2_global)*A2(20)*S2p_local ...
         + 2*vec_tilde(w2_global)*A2(20)*S2p_vel_local ...
         + A2(20)*[0 0 0]'
     
     
%

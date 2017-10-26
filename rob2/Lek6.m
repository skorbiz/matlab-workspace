%%
% Problem 2
%
clc;clear;close all
format compact

alpha = 50*pi/180

A = [cos(alpha) sin(alpha); -sin(alpha) cos(alpha)]

I=A*A'

P_local = [80; -55]
P_global = A'*P_local

Q_global = [-30; -60]
Q_local = A*Q_global

w = 0.5
Ad = [-sin(alpha) cos(alpha); -cos(alpha) -sin(alpha)]
v = w*Ad'*P_local


%%
% Problem 3
% 
clc;clear;close all
format compact

%function y = A(x)
%x = x*pi/180
%y = [cos(x) sin(x); -sin(x) cos(x)]
%end

d2r = pi/180;
A = @(x) [ cos(x*d2r) sin(x*d2r); -sin(x*d2r)  cos(x*d2r)];          %# An anonymous function
B = @(x) [-sin(x*d2r) cos(x*d2r); -cos(x*d2r) -sin(x*d2r)];          %# An anonymous function


disp(' ')
disp('1) Calculate p_global')

p_local = [150+100*cos(30*d2r); 0 + 100*sin(30*d2r)];

s1q_local = [250; -100];
s1q_global = A(60)'*s1q_local;

p_global = s1q_global + A(0)'*p_local


disp(' ')
disp('2) Calculate w2_global')
w1_local = 2;
w2_local = -8;
w2_global = w1_local + w2_local


disp(' ')
disp('3) Find p_velocity_global')
p_velocity_global = w1_local * B(60)'*s1q_local + w2_global*B(0)'*p_local


disp(' ')
disp('4) Find p_acceleration_global')

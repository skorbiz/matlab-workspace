%% 5.12 a

a = [-1 0 0; 0.3 -0.1 0.05; 1 0 0];
c = [0 0 -2];

sys = ss(a,[],c,[]);

rank(obsv(sys))


%% 5.12 b

a = [0 0.1 -100 4; -250 -7 3 50; 0 0 -3.3 0.06; 2 0 0 0.25];
c = [0 0 1 0];

sys = ss(a,[],c,[]);
rank(obsv(sys))

%% 5.12 c

a = [0 0; 0 -1]
c = [1 0; -2 0]

sys = ss(a, [], c, [])
rank(obsv(sys))

%% 5.17

clear
a = [0.4158 1.025 -0.00267 -0.0001106 -0.08021 0; -5.5 -0.8302 -0.06549 0.0039 -5.115 0.809; 0 0 0 1.0 0 0; -1040 -78.35 -34.83 -0.6214 -865.6 -631; 0 0 0 0 -75 0; 0 0 0 0 0 -100];
b = [0 0; 0 0; 0 0; 0 0; 75 0; 0 100];
c = [-1491 -146.43 -40.2 -0.9241 -1285 -564.66; 0 1.0 0 0 0 0];
d = [0 0; 0 0];

sys = ss(a,b,c,d);
%% a) Is the system oberserverble with both inputs

rank(obsv(sys))

%% b) Is the system observerble with only y1

sysOneInput = ss(a,b,[-1491 -146.43 -40.2 -0.9241 -1285 -564.66],[]);
rank(obsv(sysOneInput))

%% c) Is the system observerble with only y2

sysOneInput = ss(a,b,[0 1.0 0 0 0 0],[]);
rank(obsv(sysOneInput))

% %% d) Design a full observerble input

p1 = -4 +4i;
p2 = -4 -4i;
p3 = -3 +3i;
p4 = -3 -3i;
p5 = -100;
p6 = -75;
p = [p1 p2 p3 p4 p5 p6];

k = place(a, b, p);

ACL = a-b*k;    %Den nye A

eig(ACL)
sysOneInput = ss(ACL,b,[-1491 -146.43 -40.2 -0.9241 -1285 -564.66],[]);
rank(obsv(sysOneInput))



% %% e) 

x0 = [0 0.5 0 0 0 0 0 0 0 0 0 0]';

p1 = -3 +3i;
p2 = -3 -3i;
p3 = -1 +2i;
p4 = -1 -2i;
p5 = -100;
p6 = -75;
p = [p1 p2 p3 p4 p5 p6];

k = place(a, b, p);

l = place(a', c(1,:)', p)';

ACL = [a-b*k b*k; zeros(6) a-l* c(1,:)]

sys = ss(ACL, zeros(12,1), eye(12), zeros(12,1));

dt = 0.01;
t = 0:dt:3;

eig(ACL)

[y,t,e]=initial(sys, x0 ,t);
plot(t,e)

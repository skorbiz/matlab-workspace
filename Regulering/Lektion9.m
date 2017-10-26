%% Opgave 1
clear
a = 1;
b = -2;
c = -2;
d = 4;

A = [a b; c d];
x0 = [100 50]';

% Bestem state transition matrix (STM)
syms s;
STMs = inv((s*eye(size(A))-A))
STMt = ilaplace(STMs)*x0

% Bestem input responset for x1 = 100 og x2 =50

B= [0 0]';

sys = ss(A,B,[],[]);
t =0:0.0001:10;
u=ones(size(t));

[y,t,x]=lsim(sys,u,t,x0);
figure(1)
plot(t,x)

% Another solution
figure(2)
%subs(STMt,0)
ezplot(STMt(1),[0 10])
hold on
ezplot(STMt(2),[0 10])
title('response')
hold off

%% Opgave 2
clear
a = -8;
b = 4;
c = 4;
d = -2;

A = [-8 4; 4 -2];
x0 = [2 1]';

syms s;
STMs = inv((s*eye(size(A))-A));
STMt = ilaplace(STMs)*x0;

B = [0 0]';

sys = ss(A,B,[],[]);
t = 0:0.0001:1;
u = ones(size(t));

[y,t,x] = lsim(sys,u,t,x0);
%plot(t,x);



eq1 = STMt(1,1)
eq2 = STMt(2,1)
%solve('6/(5*exp(10*t)) + 4/5 = 8/5 - 3/(5*exp(10*t))')


eq = eq1 + '0 = 0' + eq2
solve(eq)

%% Opgave 3
clear
A = [-8 0; 0 -2];
B = [4 1]';


p = ctrb(A, B);
rank(p)

%% Opgave 4
clear
num = [1];
den = [1 2 1];
TF = tf(num,den);

[a b c d] = tf2ss(num, den);

p1 = -7;
p2 = -7;

rank(ctrb(a,b)); %Matrixen er contralable
eig(a);          %Alle polerne er negative

k = acker(a, b, [p1 p2]);

ACL = a-b*k;    %Den nye A
eig(ACL);       %Polerne er stadig negative

%% Opgave 5
clear

A = [-1 -3 0; 0 1 1; 0 2 0];
B = [1 2/3; 0 1/3; 0 1/3];

rank(ctrb(A,B))

sys = ss(A,B,[],[]);

cjordan = canon(sys, 'modal') %Jordan representation


a = [-1 0; 0 2]
b = [0.007812 0.007812; 0 0.7454;]


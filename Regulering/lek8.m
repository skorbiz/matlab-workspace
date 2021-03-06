%% opgave 4.2
clear
a = [-1.5 0.2; 0.13 0];
b = [1 0; 1 -1];
c = [];
d = [];

sys = ss(a,b,c,d);

dt = 0.1;
t = 0:dt:5;

u1 = [t]';
u2 = ones(size(t))';
u = [u1 u2];

x0 =[10; 2];
[y, t, x] = lsim(sys, u, t, x0);
plot(t,x)

%% 4.2 teachers solution
clear
syms s;
A = [-1.5 0.2; 0.13 0];
B = [1 0; 1 -1];

%Solve exsact natural response
STM =inv(s*eye(2)-A);
STMx0 = STM*[10;2];
STMt = ilaplace(STMx0);
simplify(STMt);

%then look at forced response
U =[1/s^2; 1/s];
STMB = [STM*B*U];
STMBt = ilaplace(STMB);

%Total response to u input
Xt = STMt + STMBt;
simplify(Xt);

%calc start value
subs(Xt(1),0);
subs(Xt(2),0);

%plot some values (remember they are sumbolic)
ezplot(Xt(1), [0 5])
hold on
ezplot(Xt(2), [0 5])
title('respnse for SS, x1')
hold off

%% opgave 4.6
clear
syms t;
a = [-5/13 1; 1 -5/13];
b = [exp(t/2) exp(-t)];
c = [1 1];
d = 0;

%sys = ss(a,b,c,d);

%dt = 0.1;
%t = 0:dt:5;

u1 = [t]';
u2 = ones(size(t))';
u = [u1 u2];

%x0 =[-1; 1];
%for t = 0.1:10
x = exp(a*0.1)*x(t*10)+b*u*0.1
%end

%% opgave 5.1 a
clear
a = [0 1; 1 -2];
b = [0 0; -1 1];


p =ctrb(a, b);
rank(p)
%Controllable because rank is same order as the system.


%% opgave 5.1 b
clear
a = [0 1 0; 25 0 -1; 0 0 0];
b = [0 1 5]';

p = ctrb(a, b);
det(p)
%Not controllable because determiand is equal to 0

%% opgave 5.2
clear
% Er systemerne fra opgave 5.1 stabile  
% Opgave 5.1 a er kontrolerbar, og derfor ogs� stabil, det er derfor kun
% intteresant at se p� 5.1 b
a = [0 1 0; 25 0 -1; 0 0 0];
b = [0 1 5]';

% A' = TAT'
sys = ss(a,b,[],[]);
sys = canon(sys,'modal');

%Vi kan ikke p�virke x1, da U1 er lig 0. Hvis x1 dog i sig selv er stabil
%og henfalder, s� kan vi stadig kalde systemet stabilt.
%Vi ser derfor p� pol-placeringen, og ser det ikke er stabilt.
pzplot(sys)

%Ellers kan vi sige at et system er stabiliserbart, hvis alle egen v�rdier
%er negative (egen v�rdierne ligger p� tv�rs i en modal matrix)

%% opgave 5.3
% G�r systemet contolerbar. Det g�res ved at fjerne den s�jle og r�kke vi
% ikke kan styre. Vi har s� et systemt.

a = [-5 0; 0 0];
b = [0.3202; 5.004];

%Det kan kun g�res s� let pga at matrixen er diagonaliseret.
%Det ses s� p� kontrobiliteten

p = ctrb(a,b);
det(p)

%Determinanten er forskellig fra 0, derfor er systemet stabilt.

%% Opgave 5.8
clear
a = [0.4158 1.025 -0.00267 -0.0001106 -0.08021 0; -5.5 -0.8302 -0.06549 0.0039 -5.115 0.809; 0 0 0 1.0 0 0; -1040 -78.35 -34.83 -0.6214 -865.6 -631; 0 0 0 0 -75 0; 0 0 0 0 0 -100];
b = [0 0; 0 0; 0 0; 0 0; 75 0; 0 100];
c = [-1491 -146.43 -40.2 -0.9241 -1285 -564.66; 0 1.0 0 0 0 0];
d = [0 0; 0 0];

p1 = -3 +3i;
p2 = -3 -3i;
p3 = -1 +2i;
p4 = -1 -2i;
p5 = -100;
p6 = -75;

rank(ctrb(a,b)); %Matrixen er contralable
eig(a);          %Alle polerne er negative

k = place(a, b, [p1 p2 p3 p4 p5 p6])

ACL = a-b*k;    %Den nye A
eig(ACL);       %Polerne er stadig negative

dt = 0.1;
t = 0:dt:5;

u1 = -0.1*sin(10*t)';
u2 = sin(12*t)';
u = [u1 u2];

x0 = [0; 0.5; 0; 0; 0; 0];

sys = ss(a-b*k, b, c, d);

[y,t,X]=initial(sys, x0 ,t)
plot(t,y)

%L�res l�sning
X0 = [0; 0.5; 0; 0; 0; 0];
t = 0:0.1:5;
sysCL = ss(ACL,[],c,[]);           %Systemets naturlige response
[y,t,X] = initial(sysCL, X0, t);
plot(t,y)

H = tf([20 80 100 40],[1 15 30 0 0 0])

damp(H)

pzplot(H)
pole(H)
zero(H)


%Ramp function
t = 0:0.1:3;
y=t;


figure(1)
lsim(H,y,t);

figure(2)
w=logspace(-1,3);
bode(H,w);




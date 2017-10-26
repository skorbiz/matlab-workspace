Ra = 0.51;
K = 7.55*10^-4;
J = 2.14*10^-7;
b = 7.03*10^-8;

sys = tf(K , [J*Ra Ra*b+K^2]);

step(12*sys)

%Med belastning
%b = 5.47*10^-7;

%%


Ra = 69.03;
K = 4.47*10^-4;
J = 4.38*10^-11;
b = 7.78*10^-9;


sys = tf(K , [J*Ra Ra*b+K^2]);

step(15*sys)


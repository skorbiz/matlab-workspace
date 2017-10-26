%% Exercise 3.2b
nums ={[1 1] [1 0] 1};
dens = {[1 2 3] [1 3] [1 0 0 5]};
sys = tf(nums,dens)

cjordan =       canon(sys, 'modal')         %Jordan representation
%cobserver =     canon(sys, 'companion')    %Observer representation
%ccontroler = ss(sys)                       %Controler representation

%% Exercise 3.9

num = [0.25 0.25*0.15];
den = [0.01 0.1563 0.09559 0.041383 0.004245];
sys = tf(num,den);

%cjordan =       canon(sys, 'modal')         %Jordan representation
[a b c d] =    tf2ss(num, den)            %Controller companion
%cobserver =     canon(sys, 'companion')    %Observer representation

%% Exercise 3.11
a = [-0.045 0.036 -32 -2; -0.4 -3 -0.3 250; 0 0 0 1; 0.002 -0.04 0.001 -3.2];
b = [0 0.1; -30 0; 0 0; -10 0];
c = [0 0 1 0; 0 0 0 1];
d = [0 0; 0 0];
sys = ss(a,b,c,d);

sys = ss2tf(a,b,c,d,2);
[num, den] = ss2tf(a,b,c,d,2)
num_man = {[0 0 0 0.0002 0.0022];[0 0 0.0002 0.0022 0]};
sys = tf(num_man,den);
%[a,b,c,d] =    tf2ss(num, den)         %Controller companion

%cjordan =       canon(sys, 'modal')     %Jordan representation
%cobserver =     canon(sys, 'companion') %Observer representation
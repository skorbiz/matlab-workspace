% Ross chapter 8,I problems

clc;clear;close all;

% 8.3
disp('-------------------------------------------------------------------')
disp('problem 8.3')
disp('-------------------------------------------------------------------')

n = 64
sigma = 20

mu_0 = 50             % H0: mu = mu_0
                      % H1: mu <> mu_0

disp('MANUAL BEREGNING:')
mu_hat1 = 52.5
TS1 = (mu_hat1 - mu_0)/(sigma/sqrt(n))
p1 = 1 - (normcdf(TS1) - normcdf(-TS1))
disp('ztest BEREGNING:')
[hyp_reject_5pct_1,p1,ci_5pct_1] = ztest(mu_hat1*ones(1,n),mu_0,sigma)
disp('  '),disp('  ')

disp('MANUAL BEREGNING:')
mu_hat2 = 55
TS2 = (mu_hat2 - mu_0)/(sigma/sqrt(n))
p2_manual = 1 - (normcdf(TS2) - normcdf(-TS2))
disp('ztest BEREGNING:')
[hyp_reject_5pct_2,p2,ci_5pct_2] = ztest(mu_hat2*ones(1,n),mu_0,sigma)
disp('  '),disp('  ')

disp('MANUAL BEREGNING:')
mu_hat3 = 57.5
TS3 = (mu_hat3 - mu_0)/(sigma/sqrt(n))
p3_manual = 1 - (normcdf(TS3) - normcdf(-TS3))
disp('ztest BEREGNING:')
[hyp_reject_5pct_3,p3,ci_5pct_3] = ztest(mu_hat3*ones(1,n),mu_0,sigma)
disp('  '),disp('  ')



% 8.4
disp('-------------------------------------------------------------------')
disp('problem 8.4')
disp('-------------------------------------------------------------------')

pH = [8.18, 8.16, 8.17, 8.22, 8.19, 8.17, 8.15, 8.21, 8.16, 8.18];
n = 10
sigma = 0.02

mu_0 = 8.20           % H0: mu = mu_0
                      % H1: mu <> mu_0

disp('MANUAL BEREGNING:')
mu_hat = mean(pH)
TS = (mu_hat - mu_0)/(sigma/sqrt(n))

alfa1 = 0.1
z1 = norminv(1-alfa1/2,0,1)
reject_H0_alfa1 = abs(TS) > z1

disp('ztest BEREGNING:')
[reject_H0_alfa1,p1,ci_alfa1] = ztest(pH,mu_0,sigma,alfa1)
disp('  '),disp('  ')

disp('MANUAL BEREGNING:')
mu_hat = mean(pH)
TS = (mu_hat - mu_0)/(sigma/sqrt(n))

alfa2 = 0.05
z2 = norminv(1-alfa2/2,0,1)
reject_H0_alfa2 = abs(TS) > z2

disp('ztest BEREGNING:')
[reject_H0_alfa2,p2,ci_alfa2] = ztest(pH,mu_0,sigma,alfa2)
disp('  '),disp('  ')

disp('MANUAL BEREGNING:')
mu_hat = mean(pH)
TS = (mu_hat - mu_0)/(sigma/sqrt(n))

alfa3 = 0.01
z3 = norminv(1-alfa3/2,0,1)
reject_H0_alfa3 = abs(TS) > z3

disp('ztest BEREGNING:')
[reject_H0_alfa3,p3,ci_alfa3] = ztest(pH,mu_0,sigma,alfa3)
disp('  '),disp('  ')



% 8.5
disp('-------------------------------------------------------------------')
disp('problem 8.5')
disp('-------------------------------------------------------------------')

pressure = [210, 195, 197.4, 199, 198, 202, 196, 195.5];
n = 8
sigma = 5

mu_0_min = 200        % H0: mu >= mu_0_min
                      % H1: mu < mu_0_min

disp('MANUAL BEREGNING:')
mu_hat = mean(pressure)
TS = (mu_hat - mu_0_min)/(sigma/sqrt(n))

alfa1 = 0.05
z1 = norminv(1-alfa1,0,1)
reject_H0_alfa1 = TS < -z1

disp('ztest BEREGNING:')  % ztest 'left' means critical region (rejection of H0) is 'left'
[reject_H0_alfa1,p1,ci_alfa1] = ztest(pressure,mu_0_min,sigma,alfa1,'left')  
disp('  '),disp('  ')

disp('MANUAL BEREGNING:')
mu_hat = mean(pressure)
TS = (mu_hat - mu_0_min)/(sigma/sqrt(n))

alfa2 = 0.1
z2 = norminv(1-alfa2,0,1)
reject_H0_alfa1 = TS < -z2

disp('ztest BEREGNING:')  % ztest 'left' means critical region (rejection of H0) is 'left'
[reject_H0_alfa2,p2,ci_alfa2] = ztest(pressure,mu_0_min,sigma,alfa2,'left')  
disp('  '),disp('  ')



% 8.11
disp('-------------------------------------------------------------------')
disp('problem 8.11')
disp('-------------------------------------------------------------------')

n = 20
mu_hat = 105

mu_0_max = 100        % H0: mu <= mu_0_max
                      % H1: mu > mu_0_max
                      
disp('MANUAL BEREGNING:')
sigma1 = 5
TS1 = (mu_hat - mu_0_max)/(sigma1/sqrt(n))
p1 = 1 - normcdf(TS1)

disp('ztest BEREGNING:')  % ztest 'right' means critical region (rejection of H0) is 'right'
alfa_dummy = 0.05
[reject_H0_sigma1_5pct,p1,ci_sigma1_5pct] = ztest(mu_hat*ones(1,n),mu_0_max,sigma1,alfa_dummy,'right')  
disp('  '),disp('  ')

disp('MANUAL BEREGNING:')
sigma2 = 10
TS2 = (mu_hat - mu_0_max)/(sigma2/sqrt(n))
p2 = 1 - normcdf(TS2)

disp('ztest BEREGNING:')  % ztest 'right' means critical region (rejection of H0) is 'right'
alfa_dummy = 0.05
[reject_H0_sigma2_5pct,p2,ci_sigma2_5pct] = ztest(mu_hat*ones(1,n),mu_0_max,sigma2,alfa_dummy,'right')  
disp('  '),disp('  ')

disp('MANUAL BEREGNING:')
sigma3 = 15
TS3 = (mu_hat - mu_0_max)/(sigma3/sqrt(n))
p3 = 1 - normcdf(TS3)

disp('ztest BEREGNING:')  % ztest 'right' means critical region (rejection of H0) is 'right'
alfa_dummy = 0.05
[reject_H0_sigma3_5pct,p3,ci_sigma3_5pct] = ztest(mu_hat*ones(1,n),mu_0_max,sigma3,alfa_dummy,'right')  
disp('  '),disp('  ')



% 8.12
disp('-------------------------------------------------------------------')
disp('problem 8.12')
disp('-------------------------------------------------------------------')

n = 2500
sigma = 1

mu_0_min = 3          % H0: mu >= mu_0_min
                      % H1: mu < mu_0_min

disp('MANUAL BEREGNING:')
mu_hat = 2.95
TS = (mu_hat - mu_0_min)/(sigma/sqrt(n))

alfa = 0.05
z = norminv(1-alfa,0,1)
reject_H0_alfa = TS < -z

disp('ztest BEREGNING:')  % ztest 'left' means critical region (rejection of H0) is 'left'
[reject_H0_alfa,p,ci_alfa] = ztest(mu_hat*ones(1,n),mu_0_min,sigma,alfa,'left')  
disp('  '),disp('  ')



% 8.14
disp('-------------------------------------------------------------------')
disp('problem 8.14')
disp('-------------------------------------------------------------------')

n = 36
sigma_hat = 3.1

mu_0 = 24             % H0: mu = mu_0
                      % H1: mu <> mu_0

disp('MANUAL BEREGNING:')
mu_hat = 22.5
TS = (mu_hat - mu_0)/(sigma_hat/sqrt(n))
alfa = 0.05
df = n-1
t = tinv(1-alfa/2,df)
reject_H0_alfa = abs(TS) > t
disp('ttest BEREGNING:')
[h,p,ci,stats] = ttest([(mu_hat-sqrt((n-1)/n)*sigma_hat)*ones(1,n/2) (mu_hat+sqrt((n-1)/n)*sigma_hat)*ones(1,n/2)],mu_0,alfa)  % kunstige data skabes med mu_hat,sigma_hat
disp('  '),disp('  ')



% 8.15
disp('-------------------------------------------------------------------')
disp('problem 8.15')
disp('-------------------------------------------------------------------')

n = 28
sigma_hat = 0.3

mu_0 = 0.8            % H0: mu = mu_0
                      % H1: mu <> mu_0

disp('MANUAL BEREGNING:')
mu_hat = 1.0
TS = (mu_hat - mu_0)/(sigma_hat/sqrt(n))
alfa = 0.05
df = n-1
t = tinv(1-alfa/2,df)
reject_H0_alfa = abs(TS) > t
disp('ttest BEREGNING:')
[h,p,ci,stats] = ttest([(mu_hat-sqrt((n-1)/n)*sigma_hat)*ones(1,n/2) (mu_hat+sqrt((n-1)/n)*sigma_hat)*ones(1,n/2)],mu_0,alfa)  % kunstige data skabes med mu_hat,sigma_hat
disp('  '),disp('  ')



% 8.17
disp('-------------------------------------------------------------------')
disp('problem 8.17')
disp('-------------------------------------------------------------------')

n = 100
sigma_hat = 1.1

mu_0_max = 98.6       % H0: mu <= mu_0_max
                      % H1: mu > mu_0_max

disp('MANUAL BEREGNING:')
mu_hat = 98.74
TS = (mu_hat - mu_0_max)/(sigma_hat/sqrt(n))
alfa1 = 0.05
df = n-1
t1 = tinv(1-alfa1,df)
reject_H0_alfa1 = TS > t1
disp('ttest BEREGNING:')
[h,p,ci,stats] = ttest([(mu_hat-sqrt((n-1)/n)*sigma_hat)*ones(1,n/2) (mu_hat+sqrt((n-1)/n)*sigma_hat)*ones(1,n/2)],mu_0_max,alfa1,'right')  % kunstige data skabes med mu_hat,sigma_hat
disp('  '),disp('  ')

disp('MANUAL BEREGNING:')
mu_hat = 98.74
TS = (mu_hat - mu_0_max)/(sigma_hat/sqrt(n))
alfa2 = 0.01
df = n-1
t2 = tinv(1-alfa2,df)
reject_H0_alfa2 = TS > t2
disp('ttest BEREGNING:')
[h,p,ci,stats] = ttest([(mu_hat-sqrt((n-1)/n)*sigma_hat)*ones(1,n/2) (mu_hat+sqrt((n-1)/n)*sigma_hat)*ones(1,n/2)],mu_0_max,alfa2,'right')  % kunstige data skabes med mu_hat,sigma_hat
disp('  '),disp('  ')



% 8.21
disp('-------------------------------------------------------------------')
disp('problem 8.21')
disp('-------------------------------------------------------------------')

lifetimes = [237, 242, 232, 242, 248, 230, 244, 243, 254, 262, 234, 220, 225, 236, 232, 218, 228, 240];
n = 18

mu_0_min = 240        % H0: mu >= mu_0_min
                      % H1: mu < mu_0_min

disp('MANUAL BEREGNING:')
mu_hat = mean(lifetimes)
sigma_hat = std(lifetimes)
TS = (mu_hat - mu_0_min)/(sigma_hat/sqrt(n))
df = n-1
p = tcdf(TS,df)
disp('ttest BEREGNING:')
alfa = 0.05
[h,p,ci,stats] = ttest(lifetimes,mu_0_min,alfa,'left')  
disp('  '),disp('  ')



% 8.26
disp('-------------------------------------------------------------------')
disp('problem 8.26')
disp('-------------------------------------------------------------------')

voltages = [96, 98, 105, 92, 111, 114, 99, 103, 95, 101, 106, 97];
n = 12

mu_0_min = 100        % H0: mu >= mu_0_min
                      % H1: mu < mu_0_min

disp('MANUAL BEREGNING:')
mu_hat = mean(voltages)
sigma_hat = std(voltages)
TS = (mu_hat - mu_0_min)/(sigma_hat/sqrt(n))
df = n-1
p = tcdf(TS,df)
disp('ttest BEREGNING:')
alfa = 0.05
[h,p,ci,stats] = ttest(voltages,mu_0_min,alfa,'left')  
disp('  '),disp('  ')

mu_0_max = 100        % H0: mu <= mu_0_max
                      % H1: mu > mu_0_max

disp('MANUAL BEREGNING:')
mu_hat = mean(voltages)
sigma_hat = std(voltages)
TS = (mu_hat - mu_0_max)/(sigma_hat/sqrt(n))
df = n-1
p = 1 - tcdf(TS,df)
disp('ttest BEREGNING:')
alfa = 0.05
[h,p,ci,stats] = ttest(voltages,mu_0_max,alfa,'right')  
disp('  '),disp('  ')

%% Problem 4.1
%

clc;clear;close all;
format short
load('/Users/Aq/Downloads/MUST/Lek 4/Dataset/dataset_problem_4_1');

%Data extraction
X1 = X_female;
X2 = X_male;

n1 = 24;
n2 = 24;
p = 3;

%Descriptive statistics ?hat , ?hat
X1mean = mean(X1);      %Mu^hat
S1  = cov(X1);          %Sigma^hat

X2mean = mean(X2);      %Mu^hat
S2  = cov(X2);          %Sigma^hat

Sp = ((n1-1)*S1+(n2-1)*S2)/(n1+n2-2);


% Calculate test statistic
C_p = (2*p^2+3*p-1)/(6*(p+1));
C_n = 1/(n1-1)+1/(n2-1)-1/(n1+n2-2);
C = 1-C_p*C_n

T = (n1+n2-2)*log(det(Sp))-(n1-1)*log(det(S1))-(n2-1)*log(det(S2));
T = C * T

%Hypothesis test
alpha = 0.05;
chi2 = chi2inv(1-alpha,p*(p+1)/2)
P = 1-chi2cdf(T,p*(p+1)/2)


%Normal dist test
df = 3;
z_i = chi2rnd(df,1,200);
mahanalobis_dist1 = mahal(X1,X1);
mahanalobis_dist2 = mahal(X2,X2);

figure(3)
plotmatrix(X1,'r*')

figure(4)
plotmatrix(X2,'r*')

figure(1)
qqplot(z_i,mahanalobis_dist1)

figure(2)
qqplot(z_i,mahanalobis_dist2)



%% Problem 4.2
%

clc;clear;close all;
format short
load('/Users/Aq/Downloads/MUST/Lek 4/Dataset/dataset_problem_4_2');


%Data extraction
X1 = X_female;
X2 = X_male;

n1 = 45;
n2 = 45;
p = 2;

%Descriptive statistics ?hat , ?hat
X1mean = mean(X1);      %Mu^hat
S1  = cov(X1);          %Sigma^hat

X2mean = mean(X2);      %Mu^hat
S2  = cov(X2);          %Sigma^hat

mu_dif = X1mean-X2mean;
Sp = ((n1-1)*S1+(n2-1)*S2)/(n1+n2-2)


% Calculate test statistic to test if covariance are eqaul
C_p = (2*p^2+3*p-1)/(6*(p+1));
C_n = 1/(n1-1)+1/(n2-1)-1/(n1+n2-2);
C = 1-C_p*C_n

T = (n1+n2-2)*log(det(Sp))-(n1-1)*log(det(S1))-(n2-1)*log(det(S2));
T = C * T

%Hypothesis test to test if covariances are equal
alpha = 0.05;
chi2 = chi2inv(1-alpha,p*(p+1)/2)
P = 1-chi2cdf(T,p*(p+1)/2)


%Normal dist test
df = p;
z_i = chi2rnd(df,1,200);
mahanalobis_dist1 = mahal(X1,X1);
mahanalobis_dist2 = mahal(X2,X2);

if true
figure(3)
plotmatrix(X1,'r*')

figure(4)
plotmatrix(X2,'r*')

figure(1)
qqplot(z_i,mahanalobis_dist1)

figure(2)
qqplot(z_i,mahanalobis_dist2)
end


% Test for equality of the mean vectors for female and male populations

mu0 = [0 0];
t2 = (mu_dif-mu0)*inv((1/n1+1/n2)*Sp)*(mu_dif-mu0)'


%Hypothesis test
alpha = 0.05;

t2converted = (n1+n2-p-1)/(p*(n1+n2-2))*t2; 

P = 1-fcdf(t2converted,p,n1+n2-p-1)






%Confidense region




%%  %%  %%  %%  %%  %%  %%  %%  %%  %%  %%  %%  %%  %%  %%

mu_hat_female = X1mean;
mu_hat_male = X2mean;
Sigma_hat_pooled = Sp;
n_female = n1;
n_male = n2;

mu_diff_hat = mu_hat_female-mu_hat_male
mu_diff1_hat = mu_diff_hat(1);
mu_diff2_hat = mu_diff_hat(2);
sigma1_hat = sqrt(Sigma_hat_pooled(1,1));
sigma2_hat = sqrt(Sigma_hat_pooled(2,2));

y1 = mu_diff1_hat-sqrt(Sigma_hat_pooled(1,1)):sqrt(Sigma_hat_pooled(1,1))/50:mu_diff1_hat+sqrt(Sigma_hat_pooled(1,1)); 
y2 = mu_diff2_hat-sqrt(Sigma_hat_pooled(2,2)):sqrt(Sigma_hat_pooled(2,2))/50:mu_diff2_hat+sqrt(Sigma_hat_pooled(2,2));

[Y1,Y2] = meshgrid(y1,y2);

sizeY1 = size(Y1)
sizeY2 = size(Y2)

F = mvnpdf([Y1(:) Y2(:)],mu_diff_hat,Sigma_hat_pooled);
F = reshape(F,length(y2),length(y1));

alpha = 0.05;

%c = (p*(n-1)/(n*(n-p)))*finv(1-alpha,p,n-p);
c = (p*(n_female+n_male-2)/(n_female+n_male-p-1))*(1/n_female+1/n_male)*finv(1-alpha,p,n_female+n_male-p-1);
C = (1/(2*pi*sqrt(det(Sigma_hat_pooled))))*exp(-c/2);
C = [0.0012 0.0015 0.0017 0.0018 0.0019]
contour(y1,y2,F,C)

%contour(y1,y2,F,C,'Color','k','LineWidth',3),grid,axis equal,xlabel('{\mu}_f_e_m_a_l_e_,_1 - {\mu}_m_a_l_e_,_1','Fontsize',18),ylabel('{\mu}_f_e_m_a_l_e_,_2 - {\mu}_m_a_l_e_,_2','Fontsize',18),title({'95% Confidence region and intervals for ({\mu}_f_e_m_a_l_e_,_1-{\mu}_m_a_l_e_,_1, {\mu}_f_e_m_a_l_e_,_2-{\mu}_m_a_l_e_,_2 based on observations:';'CR (black), simult.CI (red), marg.CI (blue), Bonf.CI (green)'},'Fontsize',14)

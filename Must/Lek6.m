%% Problem 6.1
%

clc;clear;close all;
format compact
format short

load('/Users/Aq/Downloads/MUST/Lek 6/Dataset/dataset_problem_6_1');

disp('-----------------------------------------')
disp('Pre calculated parameters')
disp('-----------------------------------------')

[n,r] = size(Z)

Z = [ones(n,1) Z];



disp(' ')
disp('-----------------------------------------')
disp('Unbiased Least Squares estimate')
disp('of the regression coefficients')
disp('-----------------------------------------')
beta_hat_ls = pinv(Z'*Z)*Z'*Y



disp(' ')
disp('-----------------------------------------')
disp('Unbiased estimate of Var[e]')
disp('-----------------------------------------')

H = Z*pinv(Z'*Z)*Z';

Y_hat = H*Y;

e_hat = Y-Y_hat;

SSE = sum(e_hat'*e_hat)

sigma_hat_squared = SSE/(n-(r+1))

sigma_hat = sqrt(sigma_hat_squared)



disp(' ')
disp('-----------------------------------------')
disp('SS decomposition: SST = SSR + SSE')
disp('-----------------------------------------')

Y_hat = Z*beta_hat_ls;

SST = sum((Y-mean(Y)).^2)

SSR = sum((Y_hat-mean(Y)).^2)

SSE = sum((Y-Y_hat).^2)



disp(' ')
disp('-----------------------------------------')
disp('Coefficient of determination, R_sqr,')
disp('and the adjusted coefficient of determination, R_adj_sqr')
disp('-----------------------------------------')

R_sqr = SSR/SST

R_ajd_sqr = 1 -(SSE/(n-r-1))/(SST/(n-r-1))



disp(' ')
disp('-----------------------------------------')
disp('Covariance matrix for the estimated')
disp('regression coefficients')
disp('-----------------------------------------')

sigma_hat_beta = sigma_hat_squared*pinv(Z'*Z)



disp(' ')
disp('-----------------------------------------')
disp('95% simultaneous confidence intervals')
disp('for the estimated regression coefficients')
disp('-----------------------------------------')

alpha_sim = 0.05
beta1_sim_CI = beta_hat_ls(1) + [-1 1]*sqrt((r+1)*finv(1-alpha_sim,r+1,n-r-1))*sqrt(sigma_hat_beta(1,1))
beta2_sim_CI = beta_hat_ls(2) + [-1 1]*sqrt((r+1)*finv(1-alpha_sim,r+1,n-r-1))*sqrt(sigma_hat_beta(2,2))
beta3_sim_CI = beta_hat_ls(3) + [-1 1]*sqrt((r+1)*finv(1-alpha_sim,r+1,n-r-1))*sqrt(sigma_hat_beta(3,3))

disp(' ')
disp('-----------------------------------------')
disp('95% Bonferroni confidence intervals')
disp('for the estimated regression coefficients')
disp('-----------------------------------------')

alpha_bonf = 0.05

%s = sqrt(sigma_hat_squared*pinv(Z'*Z));
%beta1_bonf_CI = beta_hat_ls(1)+[-1 1]*tinv(1-alpha_bonf/2/(r+1),n-r-1)*s(1,1)

beta1_bonf_CI = beta_hat_ls(1)+[-1 1]*tinv(1-alpha_bonf/2/(r+1),n-r-1)*sqrt(sigma_hat_beta(1,1))
beta2_bonf_CI = beta_hat_ls(2)+[-1 1]*tinv(1-alpha_bonf/2/(r+1),n-r-1)*sqrt(sigma_hat_beta(2,2))
beta3_bonf_CI = beta_hat_ls(3)+[-1 1]*tinv(1-alpha_bonf/2/(r+1),n-r-1)*sqrt(sigma_hat_beta(3,3))



disp(' ')
disp('-----------------------------------------')
disp('Extra Sum of Squares based test')
disp('H0: b1 = b2 = 0, b0 arbitrary')
disp('-----------------------------------------')

q = 0;

Z_intercept = Z(:,1:1+q);

beta_hat_intercept = pinv(Z_intercept'*Z_intercept)*Z_intercept'*Y

H_intercept = Z_intercept*pinv(Z_intercept'*Z_intercept)*Z_intercept';

Y_hat_intercept = H_intercept*Y;

SSE_intercept = sum((Y-Y_hat_intercept).^2)

T = (SSE_intercept - SSE)/(r-q)/sigma_hat_squared

alpha = 0.05;

critical_value = finv(1-alpha, r-q, n-r-1)

p_value = 1-fcdf(T, r-q, n-r-1)



disp(' ')
disp('-----------------------------------------')
disp('ANOVA table based test')
disp('-----------------------------------------')
SSR
df_SSR = r
SSE
df_SSE = n-(r+1)
test_statistic = (SSR/df_SSR)/(SSE/df_SSE)
critical_value = finv(1-alpha,df_SSR,df_SSE) 
p_value = 1-fcdf(test_statistic,df_SSR,df_SSE)



disp(' ')
disp('-----------------------------------------')
disp('A 95% prediction interval for E[Y(z0)],')
disp('the mean of new responses at z=z0=[z1 z2]=[50.5 970]')
disp('-----------------------------------------')

Z0 = [1 50.5 970];
alpha_pediction = 0.05;

Y_z0_hat = Z0*beta_hat_ls

pred_CI_mean_Y_z0 = Y_z0_hat + [-1 1]*tinv(1-alpha_pediction/2, n-r-1)*sqrt(Z0*pinv(Z'*Z)'*Z0'*sigma_hat_squared)



disp(' ')
disp('-----------------------------------------')
disp('A 95% prediction interval for Y(z0),')
disp('a single new response at z=z0')
disp('-----------------------------------------')

pred_CI_Y_z0 = Y_z0_hat




disp(' ')
disp('-----------------------------------------')
disp('Perform a model check: analysis of residuals (estimated model errors) by plotting the residuals as function of observation index, response variable and each predictor variable. Also plot a histogram and a QQ-plot for the residuals.')
disp('-----------------------------------------')

df = 2;

z_i = chi2rnd(df,1,n);

mahanalobis_dist = mahal(Z(:,2:3),Z(:,2:3));

if(true)
figure(1)
subplot(3,3,1)
plotmatrix(Z,'r*')
title('Plot of data')


subplot(3,3,2)
qqplot(z_i,mahanalobis_dist)
title('QQ-plot of Z') % OBS THERE IS ACTUALLY NO ASSUMTION OF NORMAL DISTRIBUTION HERE


subplot(3,3,4)
hist(e_hat)
title('Histogram of error')


subplot(3,3,5)
j = 1:n;
plot(j,e_hat,'r+');
line([1 n],[0 0],'color','k');
xlim([1 n]);
title('Error plot')


subplot(3,3,6)
z_i = randn(20000,1);
qqplot(e_hat,z_i)
title('QQ plot of error with assumption of a norm dist')


subplot(3,3,7)
plot(Y,e_hat,'k+');
from = min(Y);
to = max(Y);
line([from to],[0 0],'color','k');
xlim([from to]);
title('Error plot as function of Y')


subplot(3,3,8)
plot(Z(:,2),e_hat,'g+');
from = min(Z(:,2));
to = max(Z(:,2));
line([from to],[0 0],'color','k');
xlim([from to]);
title('Error plot as function of Zj1')


subplot(3,3,9)
plot(Z(:,3),e_hat,'g+');
from = min(Z(:,3));
to = max(Z(:,3));
line([from to],[0 0],'color','k');
xlim([from to]);
title('Error plot as function of Zj2')
end



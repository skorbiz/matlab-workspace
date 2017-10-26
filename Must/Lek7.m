%% Problem 7.1
%

clc;clear;close all;
format compact
format short

load('/Users/Aq/Downloads/MUST/Lek 7/Dataset/dataset_problem_7_1');

disp('-----------------------------------------')
disp('Pre calculated parameters')
disp('-----------------------------------------')

m = 4;
r = 4;
n = 62;

Z = [ones(n,1) Z];


disp(' ')
disp('-----------------------------------------')
disp('Unbiased Least Squares estimate')
disp('of the regression coefficients')
disp('-----------------------------------------')

beta_hat_ls = pinv(Z'*Z)*Z'*Y


disp(' ')
disp('-----------------------------------------')
disp('Unbiased estimate of the')
disp('(m x m) covariance matrix ?e = cov(e)')
disp('-----------------------------------------')

Y_hat = Z*beta_hat_ls;
e_hat = Y-Y_hat;
Sigma_e_hat = e_hat'*e_hat/(n-r-1)


disp(' ')
disp('-----------------------------------------')
disp('Uncorrected (not mean-corrected) decomposition')
disp('of (m x m) SS matrices: SST = SSR + SSE')
disp('-----------------------------------------')

SST_uncorrected = Y'*Y
SSE_uncorrected = e_hat'*e_hat
SSR_uncorrected = Y_hat'*Y_hat


disp(' ')
disp('-----------------------------------------')
disp('Mean-corrected decomposition')
disp('of (m x m) SS matrices: SST = SSR + SSE')
disp('-----------------------------------------')

Y_mean_matrix = ones(n,1)*mean(Y);

SST_mean_corrected = (Y-Y_mean_matrix)'*(Y-Y_mean_matrix)
SSE_mean_corrected = (Y-Y_hat)'*(Y-Y_hat)
SSR_mean_corrected = (Y_hat-Y_mean_matrix)'*(Y_hat-Y_mean_matrix)


disp(' ')
disp('-----------------------------------------')
disp('Ccoefficient of determination , R_sqr,')
disp('and the adjusted coefficient of determination, R_adj_sqr')
disp('-----------------------------------------')

R_sqr = det(SSR_mean_corrected)/det(SST_mean_corrected)

R_sqr = 1-det(SSE_mean_corrected)/det(SST_mean_corrected)

%Spørgsmål : hvorfor er de to ovenstående ikke det samme?

R_ajd_sqr = 1 -(det(SSE_mean_corrected)/(n-(r+1)*m))/(det(SST_mean_corrected)/(n-(r+1)*m))

%Spørgsmål, hvorfor pokker skald et ganges med m


disp(' ')
disp('-----------------------------------------')
disp('Inference for the regression coefficients,')
disp('indicate the ((r+1) x 1) vector of estimated regression coefficients for Y1')
disp('-----------------------------------------')

b_x1_hat = beta_hat_ls(:,1)

Sigma_hat_b_x1_hat = Sigma_e_hat(1,1)*pinv(Z'*Z)


disp(' ')
disp('-----------------------------------------')
disp('95% simultaneous confidence intervals')
disp('for the estimated regression coefficients')
disp('-----------------------------------------')

alpha_sim = 0.05
beta1_sim_CI = b_x1_hat(1) + [-1 1]*sqrt((r+1)*finv(1-alpha_sim,r+1,n-r-1))*sqrt(Sigma_hat_b_x1_hat(1,1))
beta2_sim_CI = b_x1_hat(2) + [-1 1]*sqrt((r+1)*finv(1-alpha_sim,r+1,n-r-1))*sqrt(Sigma_hat_b_x1_hat(2,2))
beta3_sim_CI = b_x1_hat(3) + [-1 1]*sqrt((r+1)*finv(1-alpha_sim,r+1,n-r-1))*sqrt(Sigma_hat_b_x1_hat(3,3))
beta4_sim_CI = b_x1_hat(4) + [-1 1]*sqrt((r+1)*finv(1-alpha_sim,r+1,n-r-1))*sqrt(Sigma_hat_b_x1_hat(4,4))
beta5_sim_CI = b_x1_hat(5) + [-1 1]*sqrt((r+1)*finv(1-alpha_sim,r+1,n-r-1))*sqrt(Sigma_hat_b_x1_hat(5,5))


disp(' ')
disp('-----------------------------------------')
disp('95% Bonferroni confidence intervals')
disp('for the estimated regression coefficients')
disp('-----------------------------------------')


alpha_bonf = 0.05
beta1_bonf_CI = beta_hat_ls(1)+[-1 1]*tinv(1-alpha_bonf/2/(r+1),n-r-1)*sqrt(Sigma_hat_b_x1_hat(1,1))
beta2_bonf_CI = beta_hat_ls(2)+[-1 1]*tinv(1-alpha_bonf/2/(r+1),n-r-1)*sqrt(Sigma_hat_b_x1_hat(2,2))
beta3_bonf_CI = beta_hat_ls(3)+[-1 1]*tinv(1-alpha_bonf/2/(r+1),n-r-1)*sqrt(Sigma_hat_b_x1_hat(3,3))
beta4_bonf_CI = beta_hat_ls(4)+[-1 1]*tinv(1-alpha_bonf/2/(r+1),n-r-1)*sqrt(Sigma_hat_b_x1_hat(4,4))
beta5_bonf_CI = beta_hat_ls(5)+[-1 1]*tinv(1-alpha_bonf/2/(r+1),n-r-1)*sqrt(Sigma_hat_b_x1_hat(5,5))






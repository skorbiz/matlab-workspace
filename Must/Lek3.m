%%
% Problem 3.1
%

%Load data
clc;clear;close all;
format short

load('/Users/Aq/Downloads/MUST/Lek 3/dataset_problem_3_1');


%Data extraction
X = obs;
size = size(X);
n = size(1,1);
p = size(1,2);


%Scatter plot 
figure(1)
plotmatrix(X,'r*')


%Descriptive statistics ?hat , ?hat
Xmean = mean(X);    %Mu^hat
Scov  = cov(X);     %Sigma^hat
Rcor  = corr(X);


%Hypothesis test
alpha = 0.05;
mu0 = [1750; 950];

t2 = (Xmean'-mu0)'*inv(Scov/n)*(Xmean'-mu0)
F  = (p*(n-1)/(n-p));
F  = F * finv(1-alpha,p,n-p)


%Calculate P-value
t2converted = (n-p)/(p*(n-1))*t2; 
pvalue_normal_approch = 1-fcdf(t2converted,p,n-p)


%Large-sample approximate test
X2 = chi2inv(1-alpha,p)
pvalue_large_sample_approach = 1-chi2cdf(t2,p)


%Confidense regions
mu_hat = Xmean';
Sigma_hat = Scov;

Sigma_hat = [0.0144 0.0117; 0.0117 0.0146]
mu_hat = [0.564; 0.603]
p = 2;
n = 42;
alpha = 0.05;

Sigma_hat_inv = inv(Sigma_hat)

[V,D] = eig(Sigma_hat)
Sigma_hat_root = V*sqrtm(D)*transpose(V)

CR = n*(mu_hat-u)' * Sigma_hat_inv * (mu_hat-u);

F = (p*(n-1))/(n-p) * finv(1-alpha, p, n-p)



%%
%Marginal confidense intervals
x1 = X(:,1);
x2 = X(:,2);

mu1_hat = mean(x1);
mu2_hat = mean(x2);

sigma1_hat = var(x1)
sigma2_hat = var(x2)

mu1_marg_ci_p = mu1_hat + tinv(1-alpha/2,n-1)*sqrt(Sigma_hat(1,1)/n);
mu1_marg_ci_m = mu1_hat - tinv(1-alpha/2,n-1)*sqrt(Sigma_hat(1,1)/n);
mu1_marg_ci = [mu1_marg_ci_m mu1_marg_ci_p]

mu2_marg_ci_p = mu2_hat + tinv(1-alpha/2,n-1)*sqrt(Sigma_hat(2,2)/n);
mu2_marg_ci_m = mu2_hat - tinv(1-alpha/2,n-1)*sqrt(Sigma_hat(2,2)/n);
mu2_marg_ci = [mu2_marg_ci_m mu2_marg_ci_p]


%Bonferroni confidense intervals
alpha_bonf = alpha / p;
mu1_bonf_ci_p = mu1_hat + tinv(1-alpha_bonf/2,n-1)*sqrt(Sigma_hat(1,1)/n);
mu1_bonf_ci_m = mu1_hat - tinv(1-alpha_bonf/2,n-1)*sqrt(Sigma_hat(1,1)/n);
mu1_bonf_ci = [mu1_bonf_ci_m mu1_bonf_ci_p]

mu2_bonf_ci_p = mu2_hat + tinv(1-alpha_bonf/2,n-1)*sqrt(Sigma_hat(2,2)/n);
mu2_bonf_ci_m = mu2_hat - tinv(1-alpha_bonf/2,n-1)*sqrt(Sigma_hat(2,2)/n);
mu2_bonf_ci = [mu2_bonf_ci_m mu2_bonf_ci_p]


%Simultanius confidense intervals
pnvar = (p*(n-1)/(n-p));

mu1_simu_ci_p = mu1_hat + finv(1-alpha/2,p,n-p)*sqrt(Sigma_hat(1,1)/n);
mu1_simu_ci_m = mu1_hat - finv(1-alpha/2,p,n-p)*sqrt(Sigma_hat(1,1)/n);
mu1_simu_ci = [mu1_simu_ci_m mu1_simu_ci_p]

%Plot of CR and CI's




%Model check to check
mahanalobis_dist = zeros(1,n);
for i = 1:n
    mahanalobis_dist(1,i) =(X(i,:)'-Xmean')' * inv(Scov) * (X(i,:)'-Xmean');
end

df = 2;
z_i = chi2rnd(df,1,n);

qqplot(mahanalobis_dist, z_i),grid,xlabel('quantiles for d_i^2','Fontsize',16),ylabel('quantiles for {\chi}_5^2 distribution','Fontsize',16),...
    title('qq-plot for d_i^2 versus {\chi}_5^2 distribution','Fontsize',16)




%% Problem 4.1
%

clc;clear;close all;
format compact
load('/Users/Aq/Downloads/MUST/Lek 5/Dataset/dataset_problem_5_1');


% Descriptive statistic
p = 4;
g = 3;

n1 = 30;
n2 = 30;
n3 = 30;
n = n1+n2+n3;

X1 = X_time1;
X2 = X_time2;
X3 = X_time3;

S1 = cov(X1);
S2 = cov(X2);
S3 = cov(X3);



% Bartlett test - Test for equal covariance.
disp('::: Barlet Test :::::::::::::::::::::')
Sp = ((n1-1)*S1+(n2-1)*S2+(n3-1)*S3)/(n1+n2+n3-3);

C =  1/(n1-1)+1/(n2-1)+1/(n3-1)-1/(n1+n2+n3-3);
C = (2*p^2+3*p-1)/(6*(p+1)*(g-1)) * C;
C =  1 -C;

LRT = (n-g)*log(det(Sp)) -(n1-1)*log(det(S1)) -(n2-1)*log(det(S2)) -(n3-1)*log(det(S3));

alpha = 0.05;

T = C*LRT
criticalValue = chi2inv(1-alpha,1/2*p*(p+1)*(g-1))

Pvalue = 1-chi2cdf(T,1/2*p*(p+1)*(g-1))

if(false)
disp('T < criticalValue')
disp('and Pvalue is above the 5%')
disp('The test is therefore accepted, as hoped')
disp('All covariances can be considered equal within the 5% confidence interval')
end



% Normal distribution test - For sample 1
disp('::: Normal distribution test ::::::::')

df = 3;
z_i = chi2rnd(df,1,200);
mahanalobis_dist1 = mahal(X1,X1);
mahanalobis_dist2 = mahal(X2,X2);
mahanalobis_dist3 = mahal(X3,X3);

if(false)
figure(1)
subplot(2,3,1)
plotmatrix(X1,'r*')
title('Population 1')

subplot(2,3,2)
plotmatrix(X2,'r*')
title('Population 2')

subplot(2,3,3)
plotmatrix(X2,'r*')
title('Population 3')

subplot(2,3,4)
qqplot(z_i,mahanalobis_dist1)

subplot(2,3,5)
qqplot(z_i,mahanalobis_dist2)

subplot(2,3,6)
qqplot(z_i,mahanalobis_dist2)
end



% Test for equality of skull size mean vectors from the three time periods (MANOVA1)
disp('::: Test for equality of skull size mean vectors from the three time periods (MANOVA1) ::::::::')

mu_hat = 1/n*(sum(X1) + sum(X2) + sum(X3));

mu_hat1 = mean(X1);
mu_hat2 = mean(X2);
mu_hat3 = mean(X3);

Tau_hat1 = mu_hat-mu_hat1;
Tau_hat2 = mu_hat-mu_hat2;
Tau_hat3 = mu_hat-mu_hat3;

%mu_hat_mat = zeros(n1,p);
%for i = 1:n1
%mu_hat_mat(i,:) = mu_hat;
%end % Same as ones(n1,1)*mu_hat1

SSt1 = (X1-ones(n1,1)*mu_hat)'*(X1-ones(n1,1)*mu_hat);
SSt2 = (X2-ones(n2,1)*mu_hat)'*(X2-ones(n2,1)*mu_hat);
SSt3 = (X3-ones(n3,1)*mu_hat)'*(X3-ones(n3,1)*mu_hat);
SSt  = SSt1+SSt2+SSt3;

SSw1 = (X1-ones(n1,1)*mu_hat1)'*(X1-ones(n1,1)*mu_hat1);
SSw2 = (X2-ones(n2,1)*mu_hat2)'*(X2-ones(n2,1)*mu_hat2);
SSw3 = (X3-ones(n3,1)*mu_hat3)'*(X3-ones(n3,1)*mu_hat3);
SSw  = SSw1+SSw2+SSw3;

SSb1 = (ones(n1,1)*mu_hat1-ones(n1,1)*mu_hat)'*(ones(n1,1)*mu_hat1-ones(n1,1)*mu_hat);
SSb2 = (ones(n2,1)*mu_hat2-ones(n2,1)*mu_hat)'*(ones(n2,1)*mu_hat2-ones(n3,1)*mu_hat);
SSb3 = (ones(n3,1)*mu_hat3-ones(n3,1)*mu_hat)'*(ones(n3,1)*mu_hat3-ones(n3,1)*mu_hat);
SSb  = SSb1+SSb2+SSb3;

Sw = SSw/(n-g);
Sb = SSb/(g-1);

Lambda = det(SSw)/(det(SSb+SSw));

%Speciel distribution applies in this manova case, as p>1 and g = 3 
disp('::: If exact distribution is used ::::::::')

alpha = 0.05;
T = ((n-p-2)/p) * ((1-sqrt(Lambda))/sqrt(Lambda))
criticalValue = finv(1-alpha,2*p,2*(n-p-2))

Pvalue = 1-fcdf(T, 2*p,2*(n-p-2))

if(false)
disp('T > criticalValue')
disp('and Pvalue is above the 5%')
disp('The test is therefore rejected, as hoped')
disp('The mean vectors cant be considerd equal within the 5% confidens interval')
end


%If only the approixmate distribution had been used
disp('::: If approximate distribution is used ::::::::')
alpha = 0.05;
T = -(n-1-(p+g)/2)*log(Lambda)
criticalValue = chi2inv(1-alpha,p*(g-1))

Pvalue = 1-chi2cdf(T, p*(g-1))
disp('Can also be rejected. Note that it says that the samples are even more unlikely to be equal')


%If the matlab builin MANOVA had been used
disp('::: If approximate distribution is used ::::::::')
alpha = 0.05;

XbuilinFormat = [X1; X2; X3];
Xgruping = [ones(n1,1); 2*ones(n2,1); 3*ones(n3,1)];

[dim,pValues,stats] = manova1(XbuilinFormat, Xgruping,alpha)

disp('Ask why two values appears')



% Confidence intervals for difference in mean vectors
uniquePairs = 1/2*g*(g-1)
variables = p
CIs = uniquePairs*variables



%Bonferroni confidense intervals
alpha_bonf_manova = alpha / (g*(g-1)*p);

disp(' ')
disp('Variable 1')
mu1_bonf_ci_12_p = (mu_hat1(1) -mu_hat2(1)) + tinv(1-alpha_bonf_manova,n-g)*sqrt(Sw(1,1)*(1/n1+1/n2));
mu1_bonf_ci_12_m = (mu_hat1(1) -mu_hat2(1)) - tinv(1-alpha_bonf_manova,n-g)*sqrt(Sw(1,1)*(1/n1+1/n2));
mu1_bonf_ci_12 = [mu1_bonf_ci_12_m mu1_bonf_ci_12_p]

mu1_bonf_ci_13_p = (mu_hat1(1) -mu_hat3(1)) + tinv(1-alpha_bonf_manova,n-g)*sqrt(Sw(1,1)*(1/n1+1/n3));
mu1_bonf_ci_13_m = (mu_hat1(1) -mu_hat3(1)) - tinv(1-alpha_bonf_manova,n-g)*sqrt(Sw(1,1)*(1/n1+1/n3));
mu1_bonf_ci_13 = [mu1_bonf_ci_13_m mu1_bonf_ci_13_p]

mu1_bonf_ci_23_p = (mu_hat2(1) -mu_hat3(1)) + tinv(1-alpha_bonf_manova,n-g)*sqrt(Sw(1,1)*(1/n2+1/n3));
mu1_bonf_ci_23_m = (mu_hat2(1) -mu_hat3(1)) - tinv(1-alpha_bonf_manova,n-g)*sqrt(Sw(1,1)*(1/n2+1/n3));
mu1_bonf_ci_23 = [mu1_bonf_ci_23_m mu1_bonf_ci_23_p]

disp(' ')
disp('Variable 2');
mu2_bonf_ci_12 = (mu_hat1(2) -mu_hat2(2)) + [-1 1]*tinv(1-alpha_bonf_manova,n-g)*sqrt(Sw(2,2)*(1/n1+1/n2))
mu2_bonf_ci_13 = (mu_hat1(2) -mu_hat3(2)) + [-1 1]*tinv(1-alpha_bonf_manova,n-g)*sqrt(Sw(2,2)*(1/n1+1/n3))
mu2_bonf_ci_23 = (mu_hat2(2) -mu_hat3(2)) + [-1 1]*tinv(1-alpha_bonf_manova,n-g)*sqrt(Sw(2,2)*(1/n2+1/n3))

disp(' ')
disp('Variable 3')
mu3_bonf_ci_12 = (mu_hat1(3) -mu_hat2(3)) + [-1 1]*tinv(1-alpha_bonf_manova,n-g)*sqrt(Sw(3,3)*(1/n1+1/n2))
mu3_bonf_ci_13 = (mu_hat1(3) -mu_hat3(3)) + [-1 1]*tinv(1-alpha_bonf_manova,n-g)*sqrt(Sw(3,3)*(1/n1+1/n3))
mu3_bonf_ci_23 = (mu_hat2(3) -mu_hat3(3)) + [-1 1]*tinv(1-alpha_bonf_manova,n-g)*sqrt(Sw(3,3)*(1/n2+1/n3))

disp(' ')
disp('Variable 4')
mu4_bonf_ci_12 = (mu_hat1(4) -mu_hat2(4)) + [-1 1]*tinv(1-alpha_bonf_manova,n-g)*sqrt(Sw(4,4)*(1/n1+1/n2))
mu4_bonf_ci_13 = (mu_hat1(4) -mu_hat3(4)) + [-1 1]*tinv(1-alpha_bonf_manova,n-g)*sqrt(Sw(4,4)*(1/n1+1/n3))
mu4_bonf_ci_23 = (mu_hat2(4) -mu_hat3(4)) + [-1 1]*tinv(1-alpha_bonf_manova,n-g)*sqrt(Sw(4,4)*(1/n2+1/n3))

disp('Skulle vi ikke være istand til at se hvor der ikke var overlap?')




%%                                              %%
%#      Lektion 2                               %%
%#                                              %%


clc;clear;close all;
format compact
load('/Users/Aq/Downloads/MUST/Lek 5/Dataset/dataset_problem_5_2');


b = 3;
l = 2;

p = 3;
n = 2;


X11 = X( 1: 2,3:5);
X12 = X( 5: 6,3:5);
X13 = X( 9:10,3:5);
X21 = X( 3: 4,3:5);
X22 = X( 7: 8,3:5);
X23 = X(11:12,3:5);


% Sample means
disp('---------------------------------------------------------------------------')
disp('Overall mean')
disp('---------------------------------------------------------------------------')
mu_hat = sum(X11+X12+X13+X21+X22+X23)/(l*b*n)


disp('---------------------------------------------------------------------------')
disp('Individual means:')
disp('---------------------------------------------------------------------------')
mu_hat_X11 = mean(X11)
mu_hat_X12 = mean(X12)
mu_hat_X13 = mean(X13)
mu_hat_X21 = mean(X21)
mu_hat_X22 = mean(X22)
mu_hat_X23 = mean(X23)


disp('---------------------------------------------------------------------------')
disp('X_hat_lx means:')
disp('---------------------------------------------------------------------------')
mu_hat_1x = 1/(b*n)*sum(X11+X12+X13)
mu_hat_2x = 1/(b*n)*sum(X21+X22+X23)


disp('---------------------------------------------------------------------------')
disp('X_hat_xk means:')
disp('---------------------------------------------------------------------------')
mu_hat_x1 = 1/(l*n)*sum(X11+X21)
mu_hat_x2 = 1/(l*n)*sum(X12+X22)
mu_hat_x3 = 1/(l*n)*sum(X13+X23)



disp('---------------------------------------------------------------------------')
disp('Calculation of Manova trick:')
disp('---------------------------------------------------------------------------')

SSw11 = (X11-ones(n,1)*mu_hat_X11)'*(X11-ones(n,1)*mu_hat_X11);
SSw12 = (X12-ones(n,1)*mu_hat_X12)'*(X12-ones(n,1)*mu_hat_X12);
SSw13 = (X13-ones(n,1)*mu_hat_X13)'*(X13-ones(n,1)*mu_hat_X13);
SSw21 = (X21-ones(n,1)*mu_hat_X21)'*(X21-ones(n,1)*mu_hat_X21);
SSw22 = (X22-ones(n,1)*mu_hat_X22)'*(X22-ones(n,1)*mu_hat_X22);
SSw23 = (X23-ones(n,1)*mu_hat_X23)'*(X23-ones(n,1)*mu_hat_X23);
SSw  = SSw11+SSw12+SSw13+SSw21+SSw22+SSw23


SSb1x = b*n*(mu_hat_1x - mu_hat)'*(mu_hat_1x - mu_hat);
SSb2x = b*n*(mu_hat_1x - mu_hat)'*(mu_hat_1x - mu_hat);
SSb_lx = SSb1x+SSb2x


SSbx1 = l*n*(mu_hat_x1 - mu_hat)'*(mu_hat_x1 - mu_hat);
SSbx2 = l*n*(mu_hat_x2 - mu_hat)'*(mu_hat_x2 - mu_hat);
SSbx3 = l*n*(mu_hat_x3 - mu_hat)'*(mu_hat_x3 - mu_hat);
SSb_xk = SSbx1+SSbx2+SSbx3


SSi11 = n*(mu_hat_X11 -mu_hat_1x -mu_hat_x1 +mu_hat)'*(mu_hat_X11-mu_hat_1x-mu_hat_x1+mu_hat);
SSi12 = n*(mu_hat_X12 -mu_hat_1x -mu_hat_x2 +mu_hat)'*(mu_hat_X12-mu_hat_1x-mu_hat_x2+mu_hat);
SSi13 = n*(mu_hat_X13 -mu_hat_1x -mu_hat_x3 +mu_hat)'*(mu_hat_X13-mu_hat_1x-mu_hat_x3+mu_hat);
SSi21 = n*(mu_hat_X21 -mu_hat_2x -mu_hat_x1 +mu_hat)'*(mu_hat_X21-mu_hat_2x-mu_hat_x1+mu_hat);
SSi22 = n*(mu_hat_X22 -mu_hat_2x -mu_hat_x2 +mu_hat)'*(mu_hat_X22-mu_hat_2x-mu_hat_x2+mu_hat);
SSi23 = n*(mu_hat_X23 -mu_hat_2x -mu_hat_x3 +mu_hat)'*(mu_hat_X23-mu_hat_2x-mu_hat_x3+mu_hat);
SSi = SSi11 + SSi12 + SSi13 + SSi21 + SSi22 + SSi23


df_W  = l*b*(n-1)
df_B1 = l-1
df_B2 = b-1
df_I  = (l-1)*(b-1)



disp('---------------------------------------------------------------------------')
disp('Test for systematic interaction and simplification of model')
disp('---------------------------------------------------------------------------')

LAMDA = det(SSw)/det(SSi + SSw)

alpha = 0.01

T = -(l*b*(n-1)-((p+1)-(l-1)*(b-1))/2)*log(LAMDA)

criticalValue = chi2inv(1-alpha, p*(l-1)*(b-1))

P_value = 1-chi2cdf(T, p*(l-1)*(b-1))

disp('Not rejected')

disp('---------------------------------------------------------------------------')
disp('Test for systematic factor 1 effect (approximate)')
disp('---------------------------------------------------------------------------')

LAMDA = det(SSw)/det(SSb_lx + SSw)

alpha = 0.01

T = -(l*b*(n-1)-((p+1)-(l-1))/2)*log(LAMDA)

criticalValue = chi2inv(1-alpha, p*(l-1))

P_value = 1-chi2cdf(T, p*(l-1))

disp('Not rejected')



disp('---------------------------------------------------------------------------')
disp('Test for systematic factor 2 effect (approximate)')
disp('---------------------------------------------------------------------------')

LAMDA = det(SSw)/det(SSb_xk + SSw)

alpha = 0.01

T = -(l*b*(n-1)-((p+1)-(b-1))/2)*log(LAMDA)

criticalValue = chi2inv(1-alpha, p*(b-1))

P_value = 1-chi2cdf(T, p*(b-1))

disp('Rejected')




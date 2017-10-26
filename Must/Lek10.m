%% Problem 10.1
%

clc;clear;close all;
format compact
format short
display_figure = false;

disp('-----------------------------------------')
disp('Loading data                             ')
disp('-----------------------------------------')

load('/Users/Aq/Downloads/MUST/Lek 10/Dataset/dataset_problem_10_1');

p = 2;
n_male = 37;
n_female = 29;

X_male = X(1:n_male,2:3);
X_female = X(n_male+1:n_male+n_female,2:3);


disp('-----------------------------------------')
disp('Descriptive statistic                    ')
disp('-----------------------------------------')

S_male = cov(X_male);
S_female = cov(X_female);

R_male = corr(X_male);
R_female = corr(X_female);

X_male_mean = mean(X_male)
X_female_mean = mean(X_female)


%display_figure = true
if(display_figure); display_figure = false;
disp('-----------------------------------------')
disp('Scatter plots                            ')
disp('-----------------------------------------')

figure()
hold on
title('Scatter Plot of male and female (blue and red)')
scatter(X_male(:,1),X_male(:,2),'b*')
scatter(X_female(:,1),X_female(:,2),'r*')
xlabel('tail length [mm]:')
ylabel('snout to vent length [mm]')
hold off

%subplot(2,3, [1 4])
%plotmatrix(X_male,'b*')
%plotmatrix(X_female,'r*')
end


%display_figure = true;
if(display_figure); display_figure = false;
disp('-----------------------------------------')
disp('Multivariate (bivar) normal model check  ')
disp('-----------------------------------------')

df = 2;

z_male_i = chi2rnd(df,1,n_male);
z_female_i = chi2rnd(df,1,n_female);

mahanalobis_dist_male = mahal(X_male,X_male);
mahanalobis_dist_female = mahal(X_female,X_female);

figure('Name','Scatter and QQ-plots of male and female data')

subplot(2,2,1)
title('Male scatter plot')
plotmatrix(X_male,'b*')

subplot(2,2,2)
title('Female scatter plot')
plotmatrix(X_female,'r*')

subplot(2,2,3)
qqplot(z_male_i,mahanalobis_dist_male)
title('QQ-plot of X_male')

subplot(2,2,4)
qqplot(z_female_i,mahanalobis_dist_female)
title('QQ-plot of X_female') 
end



disp('-----------------------------------------')
disp('Bartlett test - equal covariance matrices')
disp('-----------------------------------------')

g = 2;
n = n_male + n_female;

Sp = ((n_male-1)*S_male+(n_female-1)*S_female)/(n_male+n_female-2);

C =  1/(n_male-1)+1/(n_female-1)-1/(n_male+n_female-2);
C = (2*p^2+3*p-1)/(6*(p+1)*(g-1)) * C;
C =  1 -C;

LRT = (n-g)*log(det(Sp)) -(n_male-1)*log(det(S_male)) -(n_female-1)*log(det(S_female));

alpha = 0.05;

T = C*LRT
criticalValue = chi2inv(1-alpha,1/2*p*(p+1)*(g-1))

Pvalue = 1-chi2cdf(T,1/2*p*(p+1)*(g-1))


if(criticalValue < T)
disp(' ')
disp('T > criticalValue')
disp('and Pvalue is below the 5%')
disp('The test is therefore rejected')
disp('All covariances can  therefor NOT be considered equal within the 5% confidence interval')
else
disp(' ')
disp('T < criticalValue')
disp('and Pvalue is above the 5%')
disp('The test is therefore accepted')
disp('All covariances can therefor be considered equal within the 5% confidence interval')
end
disp(' ')


display_figure = true
if(display_figure); display_figure = false;
disp('-----------------------------------------')
disp('LDA analysis of seperation rule          ')
disp('-----------------------------------------')
disp('Note: It is tecnical incorrect to use LDA')
disp('As LDA only applyes for equal covariance ')
disp('Is done for a learning goal anyways      ')

%Assume equal cost
%Assume equal priors

[v1,v2] = meshgrid(100:220,350:700);

%LDA_lhs = @(X0) (X_male_mean-X_female_mean)*inv(Sp)*(X0-0.5*(X_male_mean+X_female_mean)')
tmp1 = (X_male_mean-X_female_mean)*inv(Sp);
tmp2 = -1/2*(X_male_mean+X_female_mean)';
LDA_lhs = @(X0) tmp1*(X0+tmp2);

[sizex,sizey] = size(v1);
d_lda = zeros(sizex, sizey);

for i = 1:351
    for j = 1:121
        x0 = [v1(i,j) v2(i,j)]';
        d_lda(i,j) = LDA_lhs(x0);
    end
end

figure()
hold on
title('Scatter Plot of male and female (blue and red)')
scatter(X_male(:,1),X_male(:,2),'b*')
scatter(X_female(:,1),X_female(:,2),'r*')
contour(v1,v2,d_lda,eps,'Color','k','LineWidth',1)
xlabel('tail length [mm]:')
ylabel('snout to vent length [mm]')
hold off
end




  %%
 %%%%                                
%%%%%%

















display_figure = true;
if(display_figure); display_figure = false;
disp('-----------------------------------------')
disp('QDA analysis of seperation rule          ')
disp('-----------------------------------------')
disp('Note: Correct as S_male != S_female      ')

[v1,v2] = meshgrid(100:220,350:700);

K = 0.5*log(det(S_male)/det(S_female)) + 0.5*(X_male_mean*inv(S_male)*X_male_mean' - X_female_mean*inv(S_female)*X_female_mean');

QDA_lhs = @(X0) -0.5*X0'*(inv(S_male) - inv(S_female))*X0 + (X_male_mean*inv(S_male) - X_female_mean*inv(S_female))*X0 - K;

%x0 = [0; 0]
%QDA_lhs(x0)

[sizex,sizey] = size(v1);
d_qda = zeros(sizex, sizey);

for i = 1:351
    for j = 1:121
        x0 = [v1(i,j) v2(i,j)]';
        d_qda(i,j) = QDA_lhs(x0);
    end
end

figure()
hold on
title('Scatter Plot of male and female (blue and red)')
scatter(X_male(:,1),X_male(:,2),'b*')
scatter(X_female(:,1),X_female(:,2),'r*')
%contour(v1,v2,d_lda,eps,'Color','y','LineWidth',1)
contour(v1,v2,d_qda,eps,'Color','k','LineWidth',1)
xlabel('tail length [mm]:')
ylabel('snout to vent length [mm]')
hold off

end


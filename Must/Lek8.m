%% Problem 8.1
%

clc;clear;close all;
format compact
format short
display_figure = true;

load('/Users/Aq/Downloads/MUST/Lek 8/Dataset/dataset_problem_8_1');


disp('-----------------------------------------')
disp('Descriptive statistics                   ')
disp('-----------------------------------------')

n = 54;
p = 8;
X_mean = mean(X_time);

disp('-----------------------------------------')
disp('Calculate descriptive statistics S and R ')
disp('-----------------------------------------')

S = cov(X_time);
R = corr(X_time);

if(display_figure)
disp('-----------------------------------------')
disp('Boxplot of the male time-data            ')
disp('-----------------------------------------')

figure(1)
subplot(2,3, [1 4])
boxplot(X_time)
title('Boxplot of times')
xlabel('Distance: 100m 200m 300m')
ylabel('Times')

disp('-----------------------------------------')
disp('Boxplot of the male time-data            ')
disp('-----------------------------------------')

subplot(2,3, [2 6])
plotmatrix(X_time,'r*')
title('Plot of data')
end

disp('-----------------------------------------')
disp('Perform a PCA on R using svd or eig      ')
disp('-----------------------------------------')

%R = [1 -2 0; -2 5 0; 0 0 2]
%p = 3

[V D] = eig(R);

[lambda_R index] = sort(diag(D),'descend');

for i = 1:p
e_R(:,i) = V(:,index(i));    
end

e_R = sign(e_R(1,1))*e_R;

e_R;
lambda_R;

%CLAUSES APPROACH
%[e_R lambda_R e_R] = svd(R);
%e_R = sign(e_R(1,1))*e_R;
%e_R
%lambda_R = diag(lambda_R)

disp('-----------------------------------------')
disp('Variance explained by PCs                ')
disp('-----------------------------------------')

PC_R_percent_explained = lambda_R/sum(lambda_R)

if(display_figure)    
figure(2)
bar(PC_R_percent_explained)
title('Pareto plot of PCs')
xlabel('PC')
ylabel('% variance explained by PCs')

for i = 1:p
    PC_percent_accumulated(i) = sum(PC_R_percent_explained(1:i));
end

hold on
plot(1:p,PC_percent_accumulated)
grid
hold off
end

disp('-----------------------------------------')
disp('Loadings for PC1 on Z and X variables    ')
disp('-----------------------------------------')

loading_Z_e1_R=e_R(:,1)'

V = diag(diag(S));
loading_X_e1_R = ((V^(-1/2))*e_R(:,1))'


disp('-----------------------------------------')
disp('Ranking of observations by PC 1          ')
disp('-----------------------------------------')

Z = zscore(X_time);

PC1_R_score_Z = Z*loading_Z_e1_R';

[time_rank_Z index_Z] = sort(PC1_R_score_Z);
time_rank_Z_10_best = index_Z(1:10)

disp('-----------------------------------------')
disp('Loadings for PC2 on Z and X variables    ')
disp('-----------------------------------------')
loading_Z_e2_R=e_R(:,2)'

V = diag(diag(S));
loading_X_e2_R = ((V^(-1/2))*e_R(:,2))'


disp('-----------------------------------------')
disp('Ranking of observations by PC 2          ')
disp('-----------------------------------------')

PC2_R_score_Z = Z*loading_Z_e2_R';

[time_rank_Z index_Z] = sort(PC2_R_score_Z);
time_rank_Z_10_best = index_Z(1:10)


disp('-----------------------------------------')
disp('Variance explained by PC 1 and 2         ')
disp('-----------------------------------------')

PC_R_percent_explained_by_PC1_and_PC2= sum(PC_R_percent_explained(1:2))


if(display_figure)    
disp('-----------------------------------------')
disp('Scatterplot of scores by PC 1 and 2      ')
disp('-----------------------------------------')

figure(3)
subplot(2,2, 1:2)
plot(PC1_R_score_Z,PC2_R_score_Z,'+')
grid
axis equal
title('ScatterPlot plot of PC1 and PC2 scores')
xlabel('PC1')
ylabel('PC2')

disp('-----------------------------------------')
disp('2Dimensional biplot for scores & loadings')
disp('-----------------------------------------')
subplot(2,2, 3)
names = char('100','200','400','800','1500','5000','10000','Marathon');
biplot(e_R(:,1:2),'scores',[PC1_R_score_Z PC2_R_score_Z],'varlabels',names)
title('2Dimensional biplot for scores & loadings')


disp('-----------------------------------------')
disp('3Dimensional biplot for scores & loadings')
disp('-----------------------------------------')
loading_Z_e3_R=e_R(:,3)'
PC3_R_score_Z = Z*loading_Z_e3_R';

subplot(2,2, 4)
names = char('100','200','400','800','1500','5000','10000','Marathon');
biplot(e_R(:,1:3),'scores',[PC1_R_score_Z PC2_R_score_Z PC3_R_score_Z],'varlabels',names)
title('3Dimensional biplot for scores & loadings')
end




%%



disp('---------------------------------------------------------------------------------')
disp('Correlation between PC 1 and 2 and original variables')
disp('---------------------------------------------------------------------------------')
Correlations_PC_1_R_time_versus_Z_or_X_variables = sqrt(LAMBDA_R_time(1,1))*PC_1_R_time_loadings_Z_variables
Correlations_PC_2_R_time_versus_Z_or_X_variables = sqrt(LAMBDA_R_time(2,2))*PC_2_R_time_loadings_Z_variables
figure('Name','Circle diagram for correlation between time-data PC 1 and 2 and original variables') 
t = 0:.01:pi;
x = cos(t);
y = sin(t);
plot(x,y,'k','LineWidth',3),axis square,title('Circle diagram for correlation between time-data PC 1 and 2 and original variables','Fontsize',14),grid,xlim([-1.2 1.2]),ylim([-1.2 1.2])
xlabel('Correlation with PC 1','Fontsize',14),ylabel('Correlation with PC 2','Fontsize',14)
hold on
plot(x,-y,'k','LineWidth',3);
plot(Correlations_PC_1_R_time_versus_Z_or_X_variables,Correlations_PC_2_R_time_versus_Z_or_X_variables,'r+','Linewidth',10)
text(0.93,0.45,'100','Fontsize',12)
text(0.75,0.38,'200','Fontsize',12)
text(0.73,0.28,'400','Fontsize',12)
text(0.77,-0.07,'800','Fontsize',12)
text(0.77,-0.15,'1500','Fontsize',12)
text(0.78,-0.24,'5000','Fontsize',12)
text(0.98,-0.27,'10000','Fontsize',12)
text(0.7,-0.31,'Marath','Fontsize',12)
text(-0.5,0.7,'Since        {\Sigma}_j_=_1_._._p r^2(X_i,PC_j) = 1','Fontsize',14)
text(-0.5,0.55,'it follows   r^2(X_i,PC_1) + r^2(X_i,PC_2) {\leq} 1','Fontsize',14)
hold off

disp('---------------------------------------------------------------------------------')
disp('Approx. simult. Bonferroni CI for variances of time-data PC 1 and 2')
disp('---------------------------------------------------------------------------------')
alpha = 0.1
psim = 2;
alpha_bonf = alpha/psim
CI_approx_lambda_1 = LAMBDA_R_time(1,1)./[1+norminv(1-alpha_bonf/2)*sqrt(2/n) 1 1-norminv(1-alpha_bonf/2)*sqrt(2/n)]
CI_approx_lambda_2 = LAMBDA_R_time(2,2)./[1+norminv(1-alpha_bonf/2)*sqrt(2/n) 1 1-norminv(1-alpha_bonf/2)*sqrt(2/n)]
disp('---------------------------------------------------------------------------------')

disp('---------------------------------------------------------------------------------')
disp('  ')
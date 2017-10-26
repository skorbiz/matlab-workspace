%% Problem 9.1
%

clc;clear;close all;
format compact
format short
display_figure = false;


load('/Users/Aq/Downloads/MUST/Lek 9/Dataset/dataset_problem_9_1');

n = 54;
p = 8;

disp('-----------------------------------------')
disp('Data Converted                           ')
disp('-----------------------------------------')

X_speed_100m = 100./(X_time(:,1));
X_speed_200m = 200./(X_time(:,2));
X_speed_400m = 400./(X_time(:,3));
X_speed_800m = 800./(X_time(:,4)*60);
X_speed_1500m = 1500./(X_time(:,5)*60);
X_speed_5000m = 5000./(X_time(:,6)*60);
X_speed_10000m = 10000./(X_time(:,7)*60);
X_speed_Marathon = 42195./(X_time(:,8)*60);

X_speed = [X_speed_100m X_speed_200m X_speed_400m...
             X_speed_800m X_speed_1500m X_speed_5000m...
             X_speed_10000m X_speed_Marathon];



disp('-----------------------------------------')
disp('Calculate descriptive statistics S and R ')
disp('-----------------------------------------')         
         
S = cov(X_speed);
R = corr(X_speed);
X_mean = mean(X_speed)


if(display_figure)
disp('-----------------------------------------')
disp('Box-plot and scatter-plot of S and R     ')
disp('-----------------------------------------')   

figure(1)
subplot(2,3, [1 4])
boxplot(X_speed)
title('Boxplot of times')
xlabel('Distance:')
ylabel('Times')

disp('-----------------------------------------')
disp('Scatter Plot of the male time-data       ')
disp('-----------------------------------------')

subplot(2,3, [2 6])
plotmatrix(X_time,'r*')
title('Plot of data')
end


disp('-----------------------------------------')
disp('PCA-based estimation of Factor model para')
disp('-----------------------------------------')

[V D] = eig(S);

[lambda_S index] = sort(diag(D),'descend');

for i = 1:p
e_S(:,i) = V(:,index(i));    
end

e_S = sign(e_S(1,1))*e_S;

e_S
lambda_S


disp('-----------------------------------------')
disp('Variance explained by PCs                ')
disp('-----------------------------------------')

PC_S_percent_explained = lambda_S/sum(lambda_S)

if(display_figure)    
figure(2)
bar(PC_S_percent_explained)
title('Pareto plot of PCs')
xlabel('PC')
ylabel('% variance explained by PCs')

for i = 1:p
    PC_percent_accumulated(i) = sum(PC_S_percent_explained(1:i));
end

hold on
plot(1:p,PC_percent_accumulated)
grid
hold off
end


disp('-----------------------------------------')
disp('Common factors needed to explain more 90%')
disp('-----------------------------------------')

m = 2


disp('-----------------------------------------')
disp('Factor model loadings in bi-plot         ')
disp('-----------------------------------------')

L_PCA = zeros(p,m);
for j = 1:m
%L_PCA(:,1) = sqrt(lambda_S(1)).*e_S(:,1);
%L_PCA(:,2) = sqrt(lambda_S(2)).*e_S(:,2);
L_PCA(:,j) = sqrt(lambda_S(j)).*e_S(:,j);
end

L_PCA
PSI_PCA = S-L_PCA*L_PCA';
PSI_PCA = diag(diag(PSI_PCA))

if(display_figure)
figure(3)
%subplot(2,2, 3)
names = char('100','200','400','800','1500','5000','10000','Marathon');
biplot(L_PCA(:,1:2),'varlabels',names)
title('2Dimensional biplot for scores & loadings')
end



disp('-----------------------------------------')
disp('Matrix of specific variances             ')
disp('-----------------------------------------')

%Maybe
%PSI_PCA

disp('-----------------------------------------')
disp('Communalities                            ')
disp('-----------------------------------------')

communality_PCA = zeros(i,1);

for i = 1:p
for j = 1:m
 
    communality_PCA(i) = communality_PCA(i) + L_PCA(i,j)*L_PCA(i,j);

end
end

communality_PCA

disp('-----------------------------------------')
disp('Variance explained by factors            ')
disp('-----------------------------------------')

variance_PCA = communality_PCA + diag(PSI_PCA)
S_PCA = diag((L_PCA'*L_PCA))


disp('-----------------------------------------')
disp('Varimax factor rotation:                 ')
disp('-----------------------------------------')

L_PCA_rotated = rotatefactors(L_PCA,'Method','varimax')

display_figure = true
if(display_figure)    

figure(4)
subplot(1,2, 1)
names = char('100','200','400','800','1500','5000','10000','Marathon');
title('2Dimensional biplot loadings')
biplot(e_S(:,1:2),'varlabels',names)


subplot(1,2, 2)

biplot(L_PCA_rotated(:,1:2),'varlabels',names)
title('2Dimensional biplot loadings (Rotated)')


end

















%%
% Problem 2.1
%
clc;clear;close all;

load('/Users/Aq/Downloads/MUST/Lek 2/dataset_problem_2_1.mat');

size = size(X);
n = size(1,1);
p = size(1,2);

Xmean = mean(X)
S = cov(X)
R = corr(X)
R = corrcoef(X)

figure(1)
plotmatrix(X,'r*')

%figure(2)
%gplotmatrix(X)

mahanalobis_dist = zeros(1,n);
for i = 1:n
    mahanalobis_dist(1,i) =(X(i,:)'-Xmean')' * inv(S) * (X(i,:)'-Xmean');
end


df = 5;
z_i = chi2rnd(df,1,200);
%qqplot(mahanalobis_dist,z_i),grid,xlabel('quantiles for d_i^2','Fontsize',16),ylabel('quantiles for {\chi}_5^2 distribution','Fontsize',16),...
%    title('qq-plot for d_i^2 versus {\chi}_5^2 distribution','Fontsize',16)

qqplot(z_i,mahanalobis_dist),grid,xlabel('quantiles for d_i^2','Fontsize',16),ylabel('quantiles for {\chi}_5^2 distribution','Fontsize',16),...
    title('qq-plot for d_i^2 versus {\chi}_5^2 distribution','Fontsize',16)

qqplot(X,X)


%%
% Problem 2.2
%
clc;clear;close all;

mu  = [1; 2; 0]
E = [4 1 2; 1 9 -3; 2 -3 5]
n = 1000;
p = 3;

[V,D] = eig(E);
Eroot = V*sqrtm(D)*transpose(V)

Y = randn(n,p);

for j=1:n,
Y(j,:) = mu + Eroot*Y(j,:)';
end

Ymean = mean(Y)
S = cov(Y)
R = corr(Y)

figure(1)
plotmatrix(Y,'r*')

figure(2)
boxplot(Y)

mahanalobis_dist = zeros(1,n);
for i = 1:n
    mahanalobis_dist(1,i) =(Y(i,:)'-Ymean')' * inv(S) * (Y(i,:)'-Ymean');
end

mahanalobis_dist = mahal(Y,Y);

figure(3) % All data combined
df = 3;
z_i = chi2rnd(df,1,n);
qqplot(mahanalobis_dist, z_i)

figure(4) % Data individualy compare to mvn-distribution
z_i = mvnrnd(mu, E, n);
qqplot(Y, z_i)


%%
% Problem 2.4
%

clc;clear;close all;

runs = 10000;
obs = zeros(runs,2);

n = 20;

a = 0;
b = 1;
c = 3;
d = 5;

u1 = (a+b)/2;
u2 = (d+c)/2;
u  = [u1 u2];

cov1 = 1/12*(b-a)^2;
cov2 = 1/12*(d-c)^2;
cov  = [cov1 0; 0 cov2];


for j=1:runs,

r1 = a + (b-a).*rand(n,1);
r2 = c + (d-c).*rand(n,1);

y = [r1 r2];

x = sqrt(n)*(mean(y)-u);

obs(j,:) = x;
end

gplotmatrix(obs)
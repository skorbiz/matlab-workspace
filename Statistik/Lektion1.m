clear
%Opgaver til Lektion 1 - 2/9-11

%% Opgave 2.16

%Data
lifetimes = [112,121,126,108,141,104,136,134,...
             121,118,143,116,108,122,127,140,...
             113,117,126,130,134,120,131,133,...
             118,125,151,147,137,140,132,119,...
             110,124,132,152,135,130,136,128];

%Mean værdi
mean(lifetimes)

%Median
median(lifetimes)

%mode
mode(lifetimes)

%Cumulativ frekves plot
cdfplot(lifetimes)

%% 2.23

dogs = [123760,48346,43575,42962,39484,36033,35388,29939,27282,22920,...
        22562,21037,20008,18218,14955,14790,14709,13312,12822,12822];

boxplot(dogs,'whisker', 10),grid

%% 2.26

grades = [3.46,3.72,3.95,3.55,3.62,3.80,3.86,3.71,3.56,3.49,...
          3.96,3.90,3.70,3.61,3.72,3.65,3.48,3.87,3.82,3.91,...
          3.69,3.67,3.72,3.66,3.79,3.75,3.93,3.74,3.50,3.83];

sample_mean = mean(grades)    %Sample mean
deviation = std(grades)     %Standart deviation

%dummy3 = find((grades > sample_mean - 1.5*deviation) & (grades < sample_mean + 1.5*deviation))
%proportion_1 = length(dummy3)/length(grades)
%proportion_1_Chebychev_lower_bound = 1 - 1/1.5^2

%% 2.30
salaries_1992 = [22340, 31825, 23153, 20108, 28902, 25040, 32603, 26596, 37951, 23145];
salaries_1993 = [22786, 32336, 23501, 20337, 29468, 25682, 33169, 27143, 39199, 23571];

scatter(salaries_1992, salaries_1993)
    ,grid,xlabel('salaries 1992'),ylabel('salaries 1993')

corrcoef(salaries_1992,salaries_1993)

%% 2.33

jan_temp = [40, 19, 41.2, 29.1, 47.8, 37.7, 48.9, 41.8, 16.1, 15.8];
july_temp = [73.2, 48.1, 81, 71.5, 62.8, 58.1, 65.7, 53.9, 58.6, 62.2];

scatter(jan_temp, july_temp)
    ,grid,xlabel('jan_temp'),ylabel('july_temp')

corrcoef(jan_temp,july_temp)

%% 2.34

%% 2.35
clear all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
filename1 = 'Data1.txt';
filename2 = 'Data2.txt';
filename3 = 'Data3.txt';
filename4 = 'Data4.txt';

%import data
data1 = readtable(filename1); 
data2 = readtable(filename2); 
data3 = readtable(filename3); 
data4 = readtable(filename4); 

d1 = table2array(data1);
d2 = table2array(data2);
d3 = table2array(data3);
d4 = table2array(data4);

% The first five moments
d1_moment = zeros(5, 1);
d2_moment = zeros(5, 1);
d3_moment = zeros(5, 1);
d4_moment = zeros(5, 1);
N = size(d1,1);
for  i = 1:5
    % Data1
    d1_moment(i, 1) = sum(d1.^i)/N;
    % Data2
    d2_moment(i, 1) = sum(d2.^i)/N;
    % Data3
    d3_moment(i, 1) = sum(d3.^i)/N;
    % Data4
    d4_moment(i, 1) = sum(d4.^i)/N;
end


% The second, third, fourth and fifth centered moments
d1_centered_moment = zeros(4, 1);
mean1 = d1_moment(1, 1);
d2_centered_moment = zeros(4, 1);
mean2 = d2_moment(1, 1);
d3_centered_moment = zeros(4, 1);
mean3 = d3_moment(1, 1);
d4_centered_moment = zeros(4, 1);
mean4 = d4_moment(1, 1);
for  i = 2:5
    % Data1
    d1_centered_moment(i, 1) = sum((d1 - mean1).^i)/N;
    % Data2
    d2_centered_moment(i, 1) = sum((d2 - mean2).^i)/N;
    % Data3
    d3_centered_moment(i, 1) = sum((d3 - mean3).^i)/N;
    % Data4
    d4_centered_moment(i, 1) = sum((d4 - mean4).^i)/N;
end

% The third (skeweness), fourth (kurtosis), fifth standardized moments 
% I can also use skeweness(d1) to calculate the third standardized moment
% and kurtosis(d1) to calculate the fourth standardized moment
% and std(d1) to calculate the deviation standard
d1_st_moment = zeros(3, 1);
std1 = sqrt(sum((d1 - mean1).^2)/N);
d2_st_moment = zeros(3, 1);
std2 = sqrt(sum((d2 - mean2).^2)/N);
d3_st_moment = zeros(3, 1);
std3 = sqrt(sum((d3 - mean3).^2)/N);
d4_st_moment = zeros(3, 1);
std4 = sqrt(sum((d4 - mean4).^2)/N);
for  i = 3:5
    % Data1
    d1_st_moment(i-2, 1) = sum(((d1 - mean1)./std1).^i)/N;
    % Data2
    d2_st_moment(i-2, 1) = sum(((d2 - mean2)./std2).^i)/N;
    % Data3
    d3_st_moment(i-2, 1) = sum(((d3 - mean3)./std3).^i)/N;
    % Data4
    d4_st_moment(i-2, 1) = sum(((d4 - mean4)./std4).^i)/N;
end

% Coefficient of Variation and Kurtosis
cv1 = std1/mean1;
cv2 = std2/mean2;
cv3 = std3/mean3;
cv4 = std4/mean4;

kurtosis1 = d1_st_moment(2, 1) - 3;
kurtosis2 = d2_st_moment(2, 1) - 3;
kurtosis3 = d3_st_moment(2, 1) - 3;
kurtosis4 = d4_st_moment(2, 1) - 3;

% The 10%, 25%, 50%, 75% and 90% percentiles
s1 = sort(d1(:, 1));
s2 = sort(d2(:, 1));
s3 = sort(d3(:, 1));
s4 = sort(d4(:, 1));

h10 = (N-1)*10/100 + 1;
h25 = (N-1)*25/100 + 1;
h50 = (N-1)*50/100 + 1;
h75 = (N-1)*75/100 + 1;
h90 = (N-1)*90/100 + 1;
h_percentiles = [h10,h25,h50,h75,h90];

d1_percentiles = zeros(5,1);
d2_percentiles = zeros(5,1);
d3_percentiles = zeros(5,1);
d4_percentiles = zeros(5,1);

for i = 1:5
    d1_percentiles(i,1) = s1(floor(h_percentiles(i)),1) + (h_percentiles(i) - floor(h_percentiles(i))) * (s1(floor(h_percentiles(i))+1,1) - s1(floor(h_percentiles(i)),1));
    d2_percentiles(i,1) = s2(floor(h_percentiles(i)),1) + (h_percentiles(i) - floor(h_percentiles(i))) * (s2(floor(h_percentiles(i))+1,1) - s2(floor(h_percentiles(i)),1));
    d3_percentiles(i,1) = s3(floor(h_percentiles(i)),1) + (h_percentiles(i) - floor(h_percentiles(i))) * (s3(floor(h_percentiles(i))+1,1) - s3(floor(h_percentiles(i)),1));
    d4_percentiles(i,1) = s4(floor(h_percentiles(i)),1) + (h_percentiles(i) - floor(h_percentiles(i))) * (s4(floor(h_percentiles(i))+1,1) - s4(floor(h_percentiles(i)),1));

end

% The cross-covariance for lags m=1, m=2 and m=3
cross_cov1 = zeros(3,1);
cross_cov2 = zeros(3,1);
cross_cov3 = zeros(3,1);
cross_cov4 = zeros(3,1);
for i = 1:3
    cross_cov1(i,1) = sum((d1(1:end-i,1)-mean1) .* (d1((i+1):end,1)-mean1))/(N-i);
    cross_cov2(i,1) = sum((d2(1:end-i,1)-mean2) .* (d2((i+1):end,1)-mean2))/(N-i);
    cross_cov3(i,1) = sum((d3(1:end-i,1)-mean3) .* (d3((i+1):end,1)-mean3))/(N-i);
    cross_cov4(i,1) = sum((d4(1:end-i,1)-mean4) .* (d4((i+1):end,1)-mean4))/(N-i);
end

% The Pearson Correlation Coefficient for lags m=1, m=2, m=3
pearson1 = zeros(3,1);
pearson2 = zeros(3,1);
pearson3 = zeros(3,1);
pearson4 = zeros(3,1);
for i = 1:3
    pearson1(i,1) = sum((d1(1:end-i,1)-mean1) .* (d1((i+1):end,1)-mean1))/(N-i)/(std1.^ 2);
    pearson2(i,1) = sum((d2(1:end-i,1)-mean2) .* (d2((i+1):end,1)-mean2))/(N-i)/(std2.^ 2);
    pearson3(i,1) = sum((d3(1:end-i,1)-mean3) .* (d3((i+1):end,1)-mean3))/(N-i)/(std3.^ 2);
    pearson4(i,1) = sum((d4(1:end-i,1)-mean4) .* (d4((i+1):end,1)-mean4))/(N-i)/(std4.^ 2);
end



disp("Data1");
disp("The first five moments:");
fprintf(1, "1st moment: %g\n", d1_moment(1,1));
fprintf(1, "2nd moment: %g\n", d1_moment(2,1));
fprintf(1, "3rd moment: %g\n", d1_moment(3,1));
fprintf(1, "4th moment: %g\n", d1_moment(4,1));
fprintf(1, "5th moment: %g\n\n", d1_moment(5,1));

disp("The centered moments:");
fprintf(1, "2nd centered moment: %g\n", d1_centered_moment(2,1));
fprintf(1, "3rd centered moment: %g\n", d1_centered_moment(3,1));
fprintf(1, "4th centeredmoment: %g\n", d1_centered_moment(4,1));
fprintf(1, "5th centered moment: %g\n\n", d1_centered_moment(5,1));

disp("The standardized moments:");
fprintf(1, "3rd standardized moment: %g\n", d1_st_moment(1,1));
fprintf(1, "4th standardized moment: %g\n", d1_st_moment(2,1));
fprintf(1, "5th standardized moment: %g\n\n", d1_st_moment(3,1));

fprintf(1, "Standard deviation: %g\n\n", std1);

fprintf(1, "Coefficient of variation: %g\n\n", cv1);

fprintf(1, "Kurtosis: %g\n\n", kurtosis1);

disp("The percentiles:");
fprintf(1, "10%% percentile: %g\n", d1_percentiles(1,1));
fprintf(1, "25%% percentile: %g\n", d1_percentiles(2,1));
fprintf(1, "50%% percentile: %g\n", d1_percentiles(3,1));
fprintf(1, "75%% percentile: %g\n", d1_percentiles(4,1));
fprintf(1, "90%% percentile: %g\n\n", d1_percentiles(5,1));

disp("The cross-covariance:");
fprintf(1, "lags m=1: %g\n", cross_cov1(1,1));
fprintf(1, "lags m=2: %g\n", cross_cov1(2,1));
fprintf(1, "lags m=3: %g\n\n", cross_cov1(3,1));

disp("The Pearson Correlation Coefficient:");
fprintf(1, "lags m=1: %g\n", pearson1(1,1));
fprintf(1, "lags m=2: %g\n", pearson1(2,1));
fprintf(1, "lags m=3: %g\n\n", pearson1(3,1));



disp("Data2");
disp("The first five moments:");
fprintf(1, "1st moment: %g\n", d2_moment(1,1));
fprintf(1, "2nd moment: %g\n", d2_moment(2,1));
fprintf(1, "3rd moment: %g\n", d2_moment(3,1));
fprintf(1, "4th moment: %g\n", d2_moment(4,1));
fprintf(1, "5th moment: %g\n\n", d2_moment(5,1));

disp("The centered moments:");
fprintf(1, "2nd centered moment: %g\n", d2_centered_moment(2,1));
fprintf(1, "3rd centered moment: %g\n", d2_centered_moment(3,1));
fprintf(1, "4th centeredmoment: %g\n", d2_centered_moment(4,1));
fprintf(1, "5th centered moment: %g\n\n", d2_centered_moment(5,1));

disp("The standardized moments:");
fprintf(1, "3rd standardized moment: %g\n", d2_st_moment(1,1));
fprintf(1, "4th standardized moment: %g\n", d2_st_moment(2,1));
fprintf(1, "5th standardized moment: %g\n\n", d2_st_moment(3,1));

fprintf(1, "Standard deviation: %g\n\n", std2);

fprintf(1, "Coefficient of variation: %g\n\n", cv2);

fprintf(1, "Kurtosis: %g\n\n", kurtosis2);

disp("The percentiles:");
fprintf(1, "10%% percentile: %g\n", d2_percentiles(1,1));
fprintf(1, "25%% percentile: %g\n", d2_percentiles(2,1));
fprintf(1, "50%% percentile: %g\n", d2_percentiles(3,1));
fprintf(1, "75%% percentile: %g\n", d2_percentiles(4,1));
fprintf(1, "90%% percentile: %g\n\n", d2_percentiles(5,1));

disp("The cross-covariance:");
fprintf(1, "lags m=1: %g\n", cross_cov2(1,1));
fprintf(1, "lags m=2: %g\n", cross_cov2(2,1));
fprintf(1, "lags m=3: %g\n\n", cross_cov2(3,1));

disp("The Pearson Correlation Coefficient:");
fprintf(1, "lags m=1: %g\n", pearson2(1,1));
fprintf(1, "lags m=2: %g\n", pearson2(2,1));
fprintf(1, "lags m=3: %g\n\n", pearson2(3,1));



disp("Data3");
disp("The first five moments:");
fprintf(1, "1st moment: %g\n", d3_moment(1,1));
fprintf(1, "2nd moment: %g\n", d3_moment(2,1));
fprintf(1, "3rd moment: %g\n", d3_moment(3,1));
fprintf(1, "4th moment: %g\n", d3_moment(4,1));
fprintf(1, "5th moment: %g\n\n", d3_moment(5,1));

disp("The centered moments:");
fprintf(1, "2nd centered moment: %g\n", d3_centered_moment(2,1));
fprintf(1, "3rd centered moment: %g\n", d3_centered_moment(3,1));
fprintf(1, "4th centeredmoment: %g\n", d3_centered_moment(4,1));
fprintf(1, "5th centered moment: %g\n\n", d3_centered_moment(5,1));

disp("The standardized moments:");
fprintf(1, "3rd standardized moment: %g\n", d3_st_moment(1,1));
fprintf(1, "4th standardized moment: %g\n", d3_st_moment(2,1));
fprintf(1, "5th standardized moment: %g\n\n", d3_st_moment(3,1));

fprintf(1, "Standard deviation: %g\n\n", std3);

fprintf(1, "Coefficient of variation: %g\n\n", cv3);

fprintf(1, "Kurtosis: %g\n\n", kurtosis3);

disp("The percentiles:");
fprintf(1, "10%% percentile: %g\n", d3_percentiles(1,1));
fprintf(1, "25%% percentile: %g\n", d3_percentiles(2,1));
fprintf(1, "50%% percentile: %g\n", d3_percentiles(3,1));
fprintf(1, "75%% percentile: %g\n", d3_percentiles(4,1));
fprintf(1, "90%% percentile: %g\n\n", d3_percentiles(5,1));

disp("The cross-covariance:");
fprintf(1, "lags m=1: %g\n", cross_cov3(1,1));
fprintf(1, "lags m=2: %g\n", cross_cov3(2,1));
fprintf(1, "lags m=3: %g\n\n", cross_cov3(3,1));

disp("The Pearson Correlation Coefficient:");
fprintf(1, "lags m=1: %g\n", pearson3(1,1));
fprintf(1, "lags m=2: %g\n", pearson3(2,1));
fprintf(1, "lags m=3: %g\n\n", pearson3(3,1));


disp("Data4");
disp("The first five moments:");
fprintf(1, "1st moment: %g\n", d4_moment(1,1));
fprintf(1, "2nd moment: %g\n", d4_moment(2,1));
fprintf(1, "3rd moment: %g\n", d4_moment(3,1));
fprintf(1, "4th moment: %g\n", d4_moment(4,1));
fprintf(1, "5th moment: %g\n\n", d4_moment(5,1));

disp("The centered moments:");
fprintf(1, "2nd centered moment: %g\n", d4_centered_moment(2,1));
fprintf(1, "3rd centered moment: %g\n", d4_centered_moment(3,1));
fprintf(1, "4th centeredmoment: %g\n", d4_centered_moment(4,1));
fprintf(1, "5th centered moment: %g\n\n", d4_centered_moment(5,1));

disp("The standardized moments :");
fprintf(1, "3rd standardized moment: %g\n", d4_st_moment(1,1));
fprintf(1, "4th standardized moment: %g\n", d4_st_moment(2,1));
fprintf(1, "5th standardized moment: %g\n\n", d4_st_moment(3,1));

fprintf(1, "Standard deviation: %g\n\n", std4);

fprintf(1, "Coefficient of variation: %g\n\n", cv4);

fprintf(1, "Kurtosis: %g\n\n", kurtosis4);

disp("The percentiles:");
fprintf(1, "10%% percentile: %g\n", d4_percentiles(1,1));
fprintf(1, "25%% percentile: %g\n", d4_percentiles(2,1));
fprintf(1, "50%% percentile: %g\n", d4_percentiles(3,1));
fprintf(1, "75%% percentile: %g\n", d4_percentiles(4,1));
fprintf(1, "90%% percentile: %g\n\n", d4_percentiles(5,1));

disp("The cross-covariance:");
fprintf(1, "lags m=1: %g\n", cross_cov4(1,1));
fprintf(1, "lags m=2: %g\n", cross_cov4(2,1));
fprintf(1, "lags m=3: %g\n\n", cross_cov4(3,1));

disp("The Pearson Correlation Coefficient:");
fprintf(1, "lags m=1: %g\n", pearson4(1,1));
fprintf(1, "lags m=2: %g\n", pearson4(2,1));
fprintf(1, "lags m=3: %g\n\n", pearson4(3,1));


% Then draw the approximated CDF of the corresponding distribution
% Compare all CDF plots
t = tiledlayout(2,2);
nexttile
plot(s1, [1:N]/N)
title('CDF Data1')
nexttile
plot(s2, [1:N]/N)
title('CDF Data2')
nexttile
plot(s3, [1:N]/N)
title('CDF Data3')
nexttile
plot(s4, [1:N]/N)
title('CDF Data4')







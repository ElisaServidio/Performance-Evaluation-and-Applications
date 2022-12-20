clear all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RESULTS FROM LOG1
filename1 = 'Log1.csv';

dA1 = readtable(filename1); %import Log1.csv
dA1 = table2array(dA1);
nA1 = size(dA1, 1);

interArr1 = sum(dA1(:,1)) / nA1; % compute Average inter-arrival time
lambda1 = 1 / interArr1; % compute arrival rate
var1 = std(dA1); % compute Variability (standard deviation of inter-arrival times)

S = 1.2; % service time 

A1 = zeros(nA1, 1); % create matrix for arrival times (to be computed)
A1(1, 1) = dA1(1); % first arrival = value of first inter-arrival time

for  i = 2:nA1
    A1(i:end, 1) = A1(i-1, 1) + dA1(i, 1); % compute arrival times
end

AC1 = [A1, zeros(nA1,1)]; % create matrix arrival and completion times

AC1(1,2) = AC1(1,1) + S; % first completion time = first arrival time + service time

for i = 1:nA1-1 % computation of completion times
    if AC1(i+1, 1) <= AC1(i,2)
       AC1(i+1, 2) = AC1(i,2) + S;   
    else
        AC1(i+1, 2) = AC1(i+1,1) + S;
    end

end

Rt1 = AC1(:,2) - AC1(:,1); 
W1 = sum(Rt1);
nC1 = nA1;
R1 = W1 /nC1; % compute Average response time

disp("LOG1 RESULTS");
fprintf(1, "Arrival rate: %g\n", lambda1);
fprintf(1, "Average inter-arrival time: %g\n", interArr1);
fprintf(1, "Variability (standard deviation of inter-arrival times): %g\n", var1);
fprintf(1, "Average response time: %g\n", R1);
% Plot correlation 1 
Correlation1 = [dA1(1:end-1,1), dA1(2:end,1)];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RESULTS FROM LOG2
filename2 = 'Log2.csv';

dA2 = readtable(filename2); %import Log2.csv
dA2 = table2array(dA2);
nA2 = size(dA2, 1);

interArr2 = sum(dA2(:,1)) / nA2;
lambda2 = 1 / interArr2;
var2 = std(dA2);

S = 1.2;

A2 = zeros(nA2, 1);
A2(1, 1) = dA2(1);

for  i = 2:nA2
    A2(i:end, 1) = A2(i-1, 1) + dA2(i, 1);
end

AC2 = [A2, zeros(nA2,1)];

AC2(1,2) = AC2(1,1) + S;

for i = 1:nA2-1
    if AC2(i+1, 1) <= AC2(i,2)
       AC2(i+1, 2) = AC2(i,2) + S;   
    else
        AC2(i+1, 2) = AC2(i+1,1) + S;
    end

end

Rt2 = AC2(:,2) - AC2(:,1);
W2 = sum(Rt2);
nC2 = nA2;
R2 = W2 /nC2;

disp("LOG2 RESULTS");
fprintf(1, "Arrival rate: %g\n", lambda2);
fprintf(1, "Average inter-arrival time: %g\n", interArr2);
fprintf(1, "Variability (standard deviation of inter-arrival times): %g\n", var2);
fprintf(1, "Average response time: %g\n", R2);
% Plot correlation 2 
Correlation2 = [dA2(1:end-1,1), dA2(2:end,1)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RESULTS FROM LOG3
filename3 = 'Log3.csv';

dA3 = readtable(filename3); %import Log3.csv
dA3 = table2array(dA3);
nA3 = size(dA3, 1);

interArr3 = sum(dA3(:,1)) / nA3;
lambda3 = 1 / interArr3;
var3 = std(dA3);

S = 1.2;

A3 = zeros(nA3, 1);
A3(1, 1) = dA3(1);

for  i = 2:nA3
    A3(i:end, 1) = A3(i-1, 1) + dA3(i, 1);
end

AC3 = [A3, zeros(nA3,1)];

AC3(1,2) = AC3(1,1) + S;

for i = 1:nA3-1
    if AC3(i+1, 1) <= AC3(i,2)
       AC3(i+1, 2) = AC3(i,2) + S;   
    else
        AC3(i+1, 2) = AC3(i+1,1) + S;
    end

end

Rt3 = AC3(:,2) - AC3(:,1);
W3 = sum(Rt3);
nC3 = nA3;
R3 = W3 /nC3;

disp("LOG3 RESULTS");
fprintf(1, "Arrival rate: %g\n", lambda3);
fprintf(1, "Average inter-arrival time: %g\n", interArr3);
fprintf(1, "Variability (standard deviation of inter-arrival times): %g\n", var3);
fprintf(1, "Average response time: %g\n", R3);
% Plot correlation 3 
Correlation3 = [dA3(1:end-1,1), dA3(2:end,1)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RESULTS FROM LOG4
filename4 = 'Log4.csv';

dA4 = readtable(filename4); %import Log4.csv
dA4 = table2array(dA4);
nA4 = size(dA4, 1);

interArr4 = sum(dA4(:,1)) / nA4;
lambda4 = 1 / interArr4;
var4 = std(dA4);

S = 1.2;

A4 = zeros(nA4, 1);
A4(1, 1) = dA4(1);

for  i = 2:nA4
    A4(i:end, 1) = A4(i-1, 1) + dA4(i, 1);
end

AC4 = [A4, zeros(nA4,1)];

AC4(1,2) = AC4(1,1) + S;

for i = 1:nA4-1
    if AC4(i+1, 1) <= AC4(i,2)
       AC4(i+1, 2) = AC4(i,2) + S;   
    else
        AC4(i+1, 2) = AC4(i+1,1) + S;
    end

end

Rt4 = AC4(:,2) - AC4(:,1);
W4 = sum(Rt4);
nC4 = nA4;
R4 = W4 /nC4;

disp("LOG4 RESULTS");
fprintf(1, "Arrival rate: %g\n", lambda4);
fprintf(1, "Average inter-arrival time: %g\n", interArr4);
fprintf(1, "Variability (standard deviation of inter-arrival times): %g\n", var4);
fprintf(1, "Average response time: %g\n", R4);
% Plot correlation 4 
Correlation4 = [dA4(1:end-1,1), dA4(2:end,1)];


% To obtain the plots visible in the delivery form I computed them individually 
% Here I display them all togheter allowing to compare them

% Compare all correlation plots
t = tiledlayout(2,2);
nexttile
plot(Correlation1(:,1), Correlation1(:,2), ".")
title('Correlation - LOG1')
nexttile
plot(Correlation2(:,1), Correlation2(:,2), ".")
title('Correlation - LOG2')
nexttile
plot(Correlation3(:,1), Correlation3(:,2), ".")
title('Correlation - LOG3')
nexttile
plot(Correlation4(:,1), Correlation4(:,2), ".")
title('Correlation - LOG4')

t.Padding = 'compact';
t.TileSpacing = 'compact';

% From comparing correlation plots:

% LOG2 - Exponential distribution
% LOG3 - Hyper-exponential distribution
% LOG4 - MMPP2 distribution

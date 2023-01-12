clear all;

%% Import and preparation of data

%import data
dataA = readtable("TraceA-A.txt"); 
dataB = readtable("TraceA-B.txt"); 
dataC = readtable("TraceA-C.txt"); 
dataD = readtable("TraceA-D.txt"); 
dataE = readtable("TraceA-E.txt"); 
A = table2array(dataA);
B = table2array(dataB);
C = table2array(dataC);
D = table2array(dataD);
E = table2array(dataE);

%preparing imported data 
N = size(A, 1);
Traces(:, 1) = A(:,1);
Traces(:, 2) = B(:,1);
Traces(:, 3) = C(:,1);
Traces(:, 4) = D(:,1);
Traces(:, 5) = E(:,1);
Traces = sort(Traces); %sort data

%% Fitting data - direct check with exponential distribution analyzing coeff of variation
%Fitting 
cv = std(Traces) ./ mean(Traces);
disp("Check if cv almost equal to 1: ")
disp(cv)
% cv almost equal to 1 --> exponential distribution
lambda = 1 ./ mean(Traces);

%compare plot of distribution and dataset
%generating exponential distribution fitting the data
T = 0:0.1:60;
for i = 1 : 5
    expCDF(i, :) = 1 - exp(- lambda(1, i) * T);
end
tl = tiledlayout(2,3);
%plot Trace A
nexttile
plot(Traces(:, 1), (1 : N) / N, T, expCDF(1, :));
title("Trace A")
legend("Dataset", "Exponential CDF",'Location','southeast')
%plot Trace B
nexttile
plot(Traces(:, 2), (1 : N) / N, T, expCDF(2, :));
title("Trace B")
legend("Dataset", "Exponential CDF",'Location','southeast')
%plot Trace C
nexttile
plot(Traces(:, 3), (1 : N) / N, T, expCDF(3, :));
title("Trace C")
legend("Dataset", "Exponential CDF",'Location','southeast')
%plot Trace D
nexttile
plot(Traces(:, 4), (1 : N) / N, T, expCDF(4, :));
title("Trace D")
legend("Dataset", "Exponential CDF",'Location','southeast')
%plot Trace E
nexttile
plot(Traces(:, 5), (1 : N) / N, T, expCDF(5, :));
title("Trace E")
legend("Dataset", "Exponential CDF",'Location','southeast')

%% Compute the minimum number of cores required to provide an average of more than 30 solutions per second

core = 0;           % multi-core machine, initialize core variable 
Tmax = 1000000;     % set maximum time of computation
solution = 0;       %solutions per second completed, at first set to 0

while solution <= 30        % compute until solution > 30
    s = 1;                  % reset initial state --> state A          
    core = core + 1;        % at each iteration add core, at first set to 1
    t = 0;                  % reset time
    sol = 0;                % solutions completed, at first set to 0
    singleCoreBusyTime = 0; % BusyTime of one core (for non parallelized stages A, C, D)
    multiCoreBusyTime = 0;  % BusyTime of all cores (for parallelized stages B, E)
    multiCoreTime = 0;      % sum of dt considering multiple cores
    while t <= Tmax
        if s == 1   %state A: create initial population
            dt = -log(rand())/lambda(1,1);
            ns = 2; %next state B
            singleCoreBusyTime = singleCoreBusyTime + dt;
        end
        if s == 2   %state B: score and scale population
            dt = (-log(rand())/lambda(1,2))/ core;
            ns = 3; %next state C
            multiCoreBusyTime = multiCoreBusyTime + dt*core;
        end
        if s == 3   %state C: retain elite
            dt = -log(rand())/lambda(1,3);
            singleCoreBusyTime = singleCoreBusyTime + dt;
            r = rand();
            if r <= 0.9
                ns = 2; %next state B
            else
                ns = 4; %next state D
            end
        end
        if s == 4   %state D: select parents
            dt = -log(rand())/lambda(1,4);
            ns = 5; %next state E
            singleCoreBusyTime = singleCoreBusyTime + dt;
        end
        if s == 5   %state E: produce crossover and mutation children
            dt = (-log(rand())/lambda(1,5))/ core;
            multiCoreBusyTime = multiCoreBusyTime + dt*core;
            r = rand();
            if r <= 0.8
                ns = 2; %next state B
            else
                ns = 1; %next state A
                sol = sol + 1; %one solution completed
            end
        end
        s = ns; %update current state
        t = t + dt;
        multiCoreTime = multiCoreTime + dt*core;
    end
    solution = sol/(t/1000); %solutions per seconds (throughput) (t from ms to s)
end

disp("Minimum number of cores required to provide an average of more than 30 solutions per second: ")
disp(core)
disp("Throughput with " + core + " core:")
disp(solution)
Utilization = (singleCoreBusyTime + multiCoreBusyTime)/multiCoreTime;
disp("Average utilization of the CPU:")
disp(Utilization)

















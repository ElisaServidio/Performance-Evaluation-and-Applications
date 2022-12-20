clear all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Consider a server, that executes jobs individually, in order of arrival and without interruption. 
% Jobs arrives and are served according to the following inter-arrival time and service time distribution:
% - Inter-arrivals: an hyper-exponential distribution with two stages, characterized by (l1 = 0.02, l2 = 0.2, p1 = 0.1)
% - Service: an uniform distribution between a=5 and b=10 

% Generate samples for the corresponding distributions and use them to create the arrival and completion curves of the service.
% Compute the 95% confidence intervals of the system average response time (R), average number of jobs (N), utilization (U) and throughput (X). 
% Consider N = 50000 jobs for the response time, and with K = 200 runs of M = 250 jobs each, for the other indices.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

NJ=50000; % Number of jobs = K*M  
K=200; % Runs
M=250; % Jobs per run

% Inter-arrival times: an hyper-exponential distribution with two stages, characterized by (l1 = 0.02, l2 = 0.2, p1 = 0.1)
inter_arrival_times = zeros(M,K);
l1 = 0.02;
l2 = 0.2;
p1 = 0.1;
% Service times: an uniform distribution between a=5 and b=10
service_times = zeros(M,K);   
a = 5;
b = 10;

% generation of services samples
service_times = a + (b-a)*(rand(M,K));

% generation of inter-arrivals samples
for j=1:K   % For each run
    for i=1:M    % For each job  
        r = rand();
        p = rand();
        if p < p1
            inter_arrival_times(i,j) = -log(r)/l1;
        else
            inter_arrival_times(i,j) = -log(r)/l2;
        end
    end
end

% Arrival times
arrival_times = zeros(M,K);
for j=1:K    % For each run
    for i=1:M    % For each job 
        if i == 1
            arrival_times(i,j) = inter_arrival_times(i,j);
        else
            arrival_times(i,j) = arrival_times(i-1,j) + inter_arrival_times(i,j);
        end
    end  
end

% Completation times
completion_times = zeros(M,K);
for j=1:K   % For each run
    for i=1:M   % For each job 
        if i == 1
            completion_times(i,j) = arrival_times(i,j) + service_times(i,j);
        else
            completion_times(i,j) = max([arrival_times(i,j), completion_times(i-1,j)]) + service_times(i,j);
        end
    end
end


response_times = completion_times - arrival_times;
% compute average response time S
R = sum(response_times)/M;
Rrow = response_times(:);

% compute average service time S
B = sum(service_times);
S = B/M;
for i = 1:K % for each run
    % compute throughput X
    X(1,i) = M/completion_times(M,i);  
    % compute utilization U
    U(1,i) = X(1,i)*S(1,i);               % utilization law
    % compute average number of jobs N
    N(1,i)= X(1,i)*R(1,i);                % little's law
end

gamma = 0.95;
d_gamma = 1.96; %percentile of the standard normal distribution corresponding to gamma = 95% 

% Compute the 95% confidence intervals of the system average response time (R)
Rm = mean(Rrow);
R2 = 1/(NJ-1)*sum((Rrow-Rm).^2);
CiR = [Rm - d_gamma *sqrt(R2/NJ), Rm + d_gamma * sqrt(R2/NJ)];

% Compute the 95% confidence intervals of the system average number of jobs (N)
Nm = mean(N);%
N2 = 1/(K-1)*sum((N-Nm).^2);
CiN = [Nm - d_gamma * sqrt(N2/K), Nm + d_gamma * sqrt(N2/K)];

% Compute the 95% confidence intervals of the system utilization (U) 
Um = mean(U);%
U2 = 1/(K-1)*sum((U-Um).^2);
CiU = [Um - d_gamma * sqrt(U2/K), Um + d_gamma * sqrt(U2/K)];

% Compute the 95% confidence intervals of the system throughput (X)
Xm = mean(X);%
X2 = 1/(K-1)*sum((X-Xm).^2);
CiX = [Xm - d_gamma * sqrt(X2/K), Xm + d_gamma * sqrt(X2/K)];

%Print values
fprintf(1, "Average system response time (R):    %f     %f\n\n", CiR(1,1), CiR(1,2));
fprintf(1, "Average system population (N):   %f      %f\n\n", CiN(1,1), CiN(1,2));
fprintf(1, "Throughput (X):  %f      %f\n\n", CiX(1,1), CiX(1,2));
fprintf(1, "Utilization (U): %f      %f\n\n", CiU(1,1), CiU(1,2));


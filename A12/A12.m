%% Performance indices of an M/M/1 and M/M/2 queue
% Consider a server that executes jobs arriving according to a Poisson process of rate l = 50 job/s, 
% and serves them with an average service time D = 15 ms. 
% Determine:
%   1. The utilization of the system
%   2. The probability of having exactly one job in the system
%   3. The probability of having less than 4 jobs in the system
%   4. The average queue length (job not in service)
%   5. The average response time
%   6. The probability that the response time is greater than 0.5 s.
%   7. The 90 percentile of the response time distribution
% After 1 year, the load has increased to l = 85 job/s, making the current solution no longer applicable. 
% The system administrator adds a second server and a load balancer: jobs enques at the load balancer, 
% and then are sent to the first available server. 
% Considering the communication time between load balancers and servers to be negligible comparted to the service times, 
% determine for this new configuration:
%   1. The average and total utilizations of the system
%   2. The probability of having exactly one job in the system
%   3. The probability of having less than 4 jobs in the system
%   4. The average queue length (job not in service)
%   5. The average response time

clear all;

% M/M/1 queue
lambda = 50; %poisson process arrival rate 
D = 15/1000; % ms --> s, average Service time (D = 1/mu)
mu = 1/D; %service rate 
rho=lambda*D; %traffic intensity = lambda/mu

disp("M/M/1 queue")
%   1. The utilization of the system
U = rho;
disp("The utilization of the system:")
disp(U)
%   2. The probability of having exactly one job in the system
p_1 = (1-rho)*(rho)^1;
disp("The probability of having exactly one job in the system:")
disp(p_1)
%   3. The probability of having less than 4 jobs in the system ( 1 - p(n>3) )
p_more3 = rho^(3+1); % The probability of having more than 3 jobs in the system
p_less4 = 1 - p_more3; % The probability of having 3 or less jobs in the system
disp("The probability of having less than 4 jobs in the system:")
disp(p_less4)
% %   4. The average queue length (job not in service)
Nq = (rho^2)/(1-rho);
disp("The average queue length (job not in service):")
disp(Nq)
%   5. The average response time
R = D/(1-rho);
disp("The average response time:")
disp(R)
%   6. The probability that the response time is greater than 0.5 s.
t = 0.5;
p_Rt = exp(-t/R);
disp("The probability that the response time is greater than 0.5 s:")
disp(p_Rt)
%   7. The 90 percentile of the response time distribution
p_90 = -log(1-90/100)*R;
disp("The 90 percentile of the response time distribution:")
disp(p_90)


% M/M/c queue - c=2
lambda = 85; %poisson process arrival rate 
rho=lambda*D/2; %traffic intensity = lambda/mu
disp("\n\nM/M/2 queue")
%   1. The average and total utilizations of the system  
U = lambda/mu;
avgU = U/2; % avgU = rho;
disp("The average and total utilizations of the system:")
disp(avgU)
disp(U)
%   2. The probability of having exactly one job in the system
p_1 = 2*(1-rho)/(1+rho)*(rho)^1;
disp("The probability of having exactly one job in the system:")
disp(p_1)
%   3. The probability of having less than 4 jobs in the system
% The probability of having exactly n job in the system
p(1) = (1-rho)/(1+rho);
for n = 1:3
    p(n+1) = 2*(1-rho)/(1+rho)*(rho)^n; % The probability of having exactly n job in the system
end
p_less4 = sum(p);
disp("The probability of having less than 4 jobs in the system:")
disp(p_less4)
%   4. The average queue length (job not in service)
% Nq = (2*rho/D)*(rho^2*D/(1-rho^2) and D = rho/lambda
Nq = 2*rho^3/(1-rho^2); 
disp("The average queue length (job not in service):")
disp(Nq)
%   5. The average response time
R = D/(1-(rho^2));
disp("The average response time:")
disp(R)
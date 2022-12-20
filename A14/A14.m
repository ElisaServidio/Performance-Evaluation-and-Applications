%% Performance indices of an G/G/c queue
% An authentication server receives jobs according to a Poisson process of rate l = 500 j/s. 
% The duration of each job is distributed according to an Hypo-Exponential, of rate µ1 = 1500 j/s and µ2 = 1000 j/s.
% Compute:
%   1. The utilization of the system
%   2. The (exact) average response time
%   3. The (exact) average number of jobs in the system
% After a year, the traffic increases and stabilizes: now it can be considered distributed according to a 4 stage Erlang distribution, with l = 4000 j/s. 
% To support this new scenario, a second authentication server is added, together with a load-balancer that holds request in a single queue, 
% and dispatches them to the first available server.
% Assuming the time required by the load balancer to be negligible (i.e., the system can be modelled with a G/G/2 queue), 
% compute:
%   1. The average utilization of the system
%   2. The approximate average response time
%   3. The approximate average number of jobs in the system

%% M/G/1 queue
disp("M/G/1 queue")
% Arrivals according to a Poisson process of rate l = 500 j/s.
lambda = 500;
% Job duration according to an Hypo-Exponential, of rate µ1 = 1500 j/s and µ2 = 1000 j/s.
mu1 = 1500;
mu2 = 1000;

% Service time
D = 1/mu1 + 1/mu2; % first order moment (mean)

rho = lambda*D;
m2 = 2*(1/mu1^2 + 1/mu2^2 + 1/(mu1*mu2)); % second order moment 

%   1. The utilization of the system
U1 = D*lambda;
disp("U1: ")
disp(U1)

%   2. The (exact) average response time
% Remaining time
w = (lambda*m2)/2;
% Queuing time
W = (rho*w)/(1-rho);
% Response time
R1 = D + W + w;
disp("R1: ")
disp(R1)

%   3. The (exact) average number of jobs in the system
N1 = lambda*R1;
disp("N1: ")
disp(N1)


%% G/G/2 queue
disp("G/G/2 queue")
% Arrivals according to a 4 stage Erlang distribution, with l = 4000 j/s
lambdaErlang = 4000;
k = 4;
T = k/lambdaErlang;
lambda = 1/T;

% Coeff. variation of the average inter-arrival time
ca = 1/sqrt(k); 

% Job duration according to an Hypo-Exponential distribution
% see results m2 and D from the M/G/1 part

% Coeff. variation of the service time
cv = sqrt(mu1^2+mu2^2)/(mu1+mu2);

rho = D/(2*T);

%   1. The average utilization of the system
U2 = rho;
disp("U2: ")
disp(U2)

%   2. The approximate average response time
R2 = D+((ca^2+cv^2)/2)*(((rho^2)*D)/(1-rho^2));
disp("R2: ")
disp(R2)

%   3. The approximate average number of jobs in the system
% Little's law
N2 = lambda*R2;
disp("N2: ")
disp(N2)
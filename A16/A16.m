%% Performance indices of a three server system
% A small system supports requests directed to both a three tier application and a file server. 
% The system components are respectively: 
%   - a web server (average service time S1 = 85 ms)
%   - a DB server (average service time S2 = 75 ms)
%   - a storage server which can also be accessed individually (average service time S3 = 60 ms). 
% All jobs arrives according to a Poisson process, with input rates 
%   - lIN[1] = 10 jobs / s to the web server,
%   - lIN[3] = 5 jobs / s to the storage server. 
% The entire system can then be modelled with the open queuing network in the figure.

% Compute:
%   1. The visits and the demands for the three stations
%   2. The utilization of the three stations
%   3. The throughput of the system
%   4. The average number of jobs in the three stations
%   5. The average system response time

%Average service times [ms]
S1 = 85;
S2 = 75;
S3 = 60;

%Input rates jobs/s
lambda_IN1 = 10;
lambda_IN3 = 5;

%   1. The visits and the demands for the three stations
% routing probabilities matrix
p11 = 0;
p12 = 1;
p13 = 0;
p21 = 0;
p22 = 0;
p23 = 1;
p31 = 0;
p32 = 0;
p33 = 0;
P = [p11,p12,p13;
     p21,p22,p23;
     p31,p32,p33];
% l[k] = lambda_IN[k]/lambda0
lambda0 = lambda_IN1 + lambda_IN3; %jobs/s, all the arrivals to the system
l = [lambda_IN1/lambda0,0,lambda_IN3/lambda0];
% visits matrix
I = eye(size(P));
vk = l/(I-P); 
disp("Vk:")
disp(vk)
% average service times
Sk = [S1,S2,S3]; 
Sk = Sk * 10^(-3); %[ms --> s]
Dk = vk .* Sk;
disp("Dk:")
disp(Dk)

%   2. The utilization of the three stations
X = lambda0; %The throughput of the system
Xk = vk * X;
Uk = Xk .* Sk;
disp("Uk:")
disp(Uk)

%   3. The throughput of the system
disp("X:")
disp(X)

%   4. The average number of jobs in the three stations
Nk = Uk./(1-Uk);
disp("Nk:")
disp(Nk)

%   5. The average system response time
Rk = Dk./(1-Uk); % Residence times 
R = sum(Rk); %average response time = sum of residence times
disp("R:")
disp(R)





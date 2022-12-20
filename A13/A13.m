%% Performance indices of an M/M/c/k system
% Consider communication channel, with a buffer of 16 packets. 
% Packets arrive at rate l = 1200 pkt/s, and require an average transmission (service) time D = 1.25 ms. 
% Modelling the system as and M/M/1/16, determine:
%   • Compute the average utilization of the channel
%   • Compute the probability of having 14 packets in the system
%   • Compute the average number of packets in the system
%   • Compute the throughput and the drop rate
%   • Compute the average response time and the average time spend in the queue
% To improve the system, and reduce the drop rate, a second channel with the same average transmission time is added. 
% The network interface sends the next packet on the first channel that becomes available. 
% Modelling the system as an M/M/2/16, compute the same performance indices:
%   • Compute the average utilization (not the total utilization) of the channel
%   • Compute the probability of having 14 packets in the system
%   • Compute the average number of packets in the system
%   • Compute the throughput and the drop rate
%   • Compute the average response time and the average time spend in the queue

%Performance indices of an M/M/1/16 system
k = 16;
lambda = 1200;
D = 1.25/1000; 
mu = 1/D;
rho = lambda/mu;

disp("Performance indices of an M/M/1/16 system")
%   • Compute the average utilization of the channel
avgU = (rho-rho^(k+1))/(1-rho^(k+1));
disp("avgU:")
disp(avgU)

%   • Compute the probability of having 14 packets in the system
p0 = (1-rho)/(1-rho^(k+1));
p_14 = p0*rho^14;
disp("p_14:")
disp(p_14)

%   • Compute the average number of packets in the system
N = rho/(1-rho)- ((k+1)*rho^(k+1))/(1-rho^(k+1));
disp("N:")
disp(N)

%   • Compute the throughput and the drop rate
p_loss = (rho^k-rho^(k+1))/(1-rho^(k+1)); 
DropRate = lambda*p_loss;
X = lambda - DropRate;
disp("X:")
disp(X)
disp("DropRate:")
disp(DropRate)

%   • Compute the average response time and the average time spend in the queue
R = N/X;
Qt = R - D;
disp("R:")
disp(R)
disp("Qt:")
disp(Qt)


%Performance indices of an M/M/2/16 system
rho = lambda/(2*mu);
c = 2;

disp("Performance indices of an M/M/2/16 system")
%   • Compute the probability of having 14 packets in the system
p(1) = ((c*rho)^c/factorial(c)*(1-rho^(k-c+1))/(1-rho)+1+(c*rho)/1)^(-1);
for i = 1: 16
    if i < c
        p(i+1) = p(1)/factorial(i)*(lambda/mu)^i;
    else
        p(i+1) = p(1)/(factorial(c)*c^(i-c))*(lambda/mu)^i;
    end 
    i = i+1;
end
disp("p14")
disp(p(15))

%   • Compute the average utilization (not the total utilization) of the channel
int = (1:c);
addend1 = sum(int.*p((int+1):(c+1)));
addend2 = c*sum(p((c+1+1):(k+1)));
U = addend1 + addend2;
avgU = U/c;
disp("avgU:")
disp(avgU)

% for i = 1:c
%    addend1 = 0;
%    addend2 = 0;
%    for j = (c+1):k
%        addend2 = addend2 + p(j+1);
%    end
%    addend1 = addend1 + i*p(i+1);
% end
% U = addend1 + c*addend2;
% avgU = U/c;

%   • Compute the average number of packets in the system
int = (0:16);
N = sum(int.*p);
disp("N");
disp(N);

%   • Compute the throughput and the drop rate
DropRate = lambda * p(k+1);
X = lambda * (1-p(k+1));
disp("X:")
disp(X)
disp("DropRate:")
disp(DropRate)

%   • Compute the average response time and the average time spent in the queue
R = N/X;
Qt = R - D;
disp("R:")
disp(R)
disp("Qt:")
disp(Qt)



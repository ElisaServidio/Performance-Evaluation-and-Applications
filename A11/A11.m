%% A communication channel with RED
% A communication protocol uses a variant of RED – Random Early Detection – to perform flow control. 

% The buffer have a size N = 16. 
% Packets arrive at rate l = 200 pkt / sec and are transmitted at a rate µ = 100 pkt / sec. 

% Whenever there are less than N0 = 8 packets in the buffer, new packets are always accepted and increase the queue length. 
% If there are more than n > N0 packets, new arrivals are discarded with probability p(n) = min(1, (n – N0) / (N – N0)) .

% Considering that all timings are exponentially distributed:
%   • Model the system with a birth death process, where the population count represents the occupation of the buffer. 
%     In particular, use the following expressions for li and µi:
%       li =    l               if i < N0
%               l*(N-i)/(N-N0)  if N0 <= i < N
%       µi = µ
%   • Determine and plot the steady state distribution of the buffer occupation [for 0 ≤ i ≤ N].
%     Please note that since the birth probability becomes zero for i ≥ N, infinite summations can be truncated at i = N


%% Probabilities of the last 7 population values p10,...,p16

N = 16;
lambda = 200; 
mu = 100;
N0 = 8;

pn = zeros(1, N+1);
pd = zeros(1, N+1);

pn(1,1) = 0;
pd(1,1) = 0;

for i = 0:N-1
    if i < N0
        lambdai = lambda;
    else
        lambdai = lambda*(N-i)/(N-N0);
    end
    mi = mu;
    pn(1, i+2) = pn(1,i+1) + log(lambdai);
    pd(1, i+2) = pd(1,i+1) + log(mi);
end 

p = exp(pn - pd);
p = p / sum(p);

disp(p(1,N-7+2:N+1))

%% Plot of the buffer length distribution
plot((0:N),p, "+", LineWidth=1);
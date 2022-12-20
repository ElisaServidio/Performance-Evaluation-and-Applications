%% Visits, demand and throughput
% An application running on an embedder systems, alternates between CPU, Disk, and network communication cycles. 
% Arrival rates, average service times, and routing probabilities have been determined and are shown in the picture.

% Compute:
% 1. The visits of the three stations
% 2. The demand of the three station
% 3. The throughput of the three stations

%% 1. The visits of the three stations
% routing probabilities
p11 = 0;
p12 = 0.3;
p13 = 0.2;
p21 = 1;
p22 = 0;
p23 = 0;
p31 = 0.2;
p32 = 0.8;
p33 = 0;
% routing probabilities matrix
P = [p11,p12,p13;
     p21,p22,p23;
     p31,p32,p33];
% λIN[k] is different from zero only for one station k
l = [1,0,0];
% visits matrix
I = eye(size(P));
v = l/(I-P); 
disp("Vk:")
disp(v)

%% 2. The demand of the three station
% average service times
S = [40,100,120]; 
S = S * 10^(-3); %[ms --> s]
D = v .* S;
disp("Dk:")
disp(D)


%% 3. The throughput of the three stations
lambda0 = 10; %jobs/s, all the arrivals to the system
X = lambda0; 
Xk = v * X;
disp("Xk:")
disp(Xk)






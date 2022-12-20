%% Communication channel availability
% A mobile application can communicate using either 4G or WiFi. 

% 4G connections is available on the average 20 hours. 
% When it goes down, it requires 2 hours to be established again. 

% Similarly, WiFi is available on the average 3 hours, and requires 8 hours to be established again. 

% To avoid software aging when both technologies are available, the mobile device resets them both every 100 hours, 
% making them both unavailable for the corresponding time previously defined.

% Assuming that all times follows the exponential distribution, 
% and that the system starts in a state where both channels are available: 

% • Draw the Markov Chain of the model
% • Compute the infinitesimal generator and solve the corresponding differential equations
% • Show the probability of the various states for the time T = [0, 300] hour

clear all;

MTTF_4G = 20;
MTTF_WiFi = 3;
MTTF_sys = 100;
MTTR_4G = 2;
MTTR_WiFi = 8;

l1 = 1/MTTF_4G;     % 1/20
l2 = 1/MTTF_WiFi;   % 1/3
l3 = 1/MTTF_sys;    % 1/100
m1 = 1/MTTR_4G;     % 1/2
m2 = 1/MTTR_WiFi;   % 1/8

Q = [-l1-l2-l3,     l2,    l1,    l3;
            m2,    -m2-l1,     0,     l1;
            m1,      0,   -m1-l2,     l2;
             0,     m1,    m2,    -m1-m2];
disp(Q(:,:))
p0 = [1, 0, 0, 0];

[t, Sol]=ode45(@(t,x) Q'*x, [0 300], p0');
plot(t, Sol, "-")
legend('s1 - 4G and WiFi up','s2 - 4G up','s3 - WiFi up','s4 - 4G and WiFi down')
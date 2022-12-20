%% A closed two tasks embedded systems
% An embedded system has four possible operational states, characterized by the following energy consumptions:
% 1. Idle [0.1 W]
% 2. CPU computation [2 W]
% 3. GPU computation [10 W]
% 4. I/O [0.5 W]

% The system starts in Idle state. 

% A New job is executed on the average every 10 sec. 

% When the job starts, it initially runs on the CPU. Then:
% • It can finish after an average time of 50 sec. and return Idle.
% • It can do I/O after an average time of 10 sec. I/O has an average duration of 5 sec, then control returns to the CPU.
% • It can do GPU computation an average time of 20 sec. GPU computation has an average duration of 2 sec, then control returns to the CPU.

% Considering that all timings are exponentially distributed, determine:
% • The probability of being in each state, both steady state and as function of time, from t =0 to t = 500;
% • The following performance metrics, both steady state and as function of time, from t =0 to t = 500:
%       o Utilization [time the system is not idle] {state reward}
%       o Average power consumption [measure in W] {state reward}
%       o System throughput [when the system returns to the idle state] {transition reward}
%       o GPU throughput [when a GPU task finishes and the control returns to the CPU] {transition reward}
%       o I/O frequency [[when a I/O task finishes and the control returns to the CPU] {transition reward}

clear all;

l1 = 1/10;     
l2 = 1/20;   
l3 = 1/10;    
m1 = 1/50;    
m2 = 1/2;   
m3 = 1/5; 

Q = [-l1,           l1,      0,      0;
      m1,    -m1-l2-l3,     l2,     l3;
       0,           m2,    -m2,      0;
       0,           m3,      0,    -m3];
p0 = [1, 0, 0, 0];

%% Probability of being in each state, both steady state and as function of time, from t =0 to t = 500
% as function of time
t =  linspace(0,500);
for i=1:size(t,2)
    pi_t(:,i) = p0 * expm(Q * t(i));
end
% steady state
u = [1, 0, 0, 0];
Qs = Q;
Qs(:,1) = ones(4,1);
pi = u * inv(Qs); 

%% State reward vectors
% 1. Idle [0.1 W]
% 2. CPU computation [2 W]
% 3. GPU computation [10 W]
% 4. I/O [0.5 W]

% Utilization [time the system is not idle] {state reward}
a1 = [0,1,1,1];
% steady state
U = sum(pi.*a1);
% as function of time
for i = 1:4
    res_U(i,:) = pi_t(i,:)*a1(i);
end 
U_t = sum(res_U);
% Average power consumption [measure in W] {state reward}
a2 = [0.1,2,10,0.5];
% steady state
power = sum(pi.*a2);
% as function of time
for i = 1:4
    res_power(i,:) = pi_t(i,:)*a2(i);
end 
power_t = sum(res_power);

%% Transition reward matrices
% System throughput [when the system returns to the idle state] {transition reward}
eps1 = [ 0,     0,      0,      0;
         1,     0,      0,      0;
         0,     0,      0,      0;
         0,     0,      0,      0];
for i = 1:4
    % steady state
    res_X(i) = pi(i)*sum(Q(i,:).*eps1(i,:));
    % as function of time
    res_X_t(i,:) = pi_t(i,:).*sum(Q(i,:).*eps1(i,:));
end
X = sum(res_X);
X_t = sum(res_X_t);
% GPU throughput [when a GPU task finishes and the control returns to the CPU] {transition reward}
eps2 = [ 0,     0,      0,      0;
         0,     0,      0,      0;
         0,     1,      0,      0;
         0,     0,      0,      0];
for i = 1:4
    % steady state
    res_X3(i) = pi(i)*sum(Q(i,:).*eps2(i,:));
    % as function of time
    res_X3_t(i,:) = pi_t(i,:).*sum(Q(i,:).*eps2(i,:));
end
X3 = sum(res_X3);
X3_t = sum(res_X3_t);
% I/O frequency [[when a I/O task finishes and the control returns to the CPU] {transition reward}
eps3 = [ 0,     0,      0,      0;
         0,     0,      0,      0;
         0,     0,      0,      0;
         0,     1,      0,      0];
for i = 1:4
    % steady state
    res_X4(i) = pi(i)*sum(Q(i,:).*eps3(i,:));
    % as function of time
    res_X4_t(i,:) = pi_t(i,:).*sum(Q(i,:).*eps3(i,:));
end
X4 = sum(res_X4); 
X4_t = sum(res_X4_t);

%% Infinitesimal generator matrix:
disp('Q:')
disp(Q)
%% State reward vectors, and transition reward matrices:
disp('a1:')
disp(a1)
disp('a2:')
disp(a2)
disp('esp1:')
disp(eps1(:,:))
disp('esp2:')
disp(eps2(:,:))
disp('esp3:')
disp(eps3(:,:))
%% Figure with the evolution of the state probabilities as function of time
tiled = tiledlayout(2,2);
nexttile
plot(t, pi_t)
legend('s1 - Idle','s2 - CPU computation','s3 - GPU computation','s4 - I/O')
%% Figure with the evolution of the rewards as function of time
nexttile
plot(t, power_t)
hold on 
plot(t, U_t)
hold off
legend('Average power consumption','Utilization')
nexttile
plot(t, X_t)
hold on 
plot(t, X3_t)
hold on 
plot(t, X4_t)
legend('System throughput', 'GPU throughput', 'I/O frequency')
%% Steady state probabilities, and limit rewards
disp('pi:')
disp(pi)
disp('power:')
disp(power)
disp('U:')
disp(U)
disp('X:')
disp(X)
disp('X3:')
disp(X3)
disp('X4:')
disp(X4)





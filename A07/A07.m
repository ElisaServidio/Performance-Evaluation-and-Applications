% An embedded system is used to control the environment condition in a room. 
% It is composed by:
%   temperature sensor                                  STATE1
%   CPU                                                 STATE2
%   actuator that can control the air conditioner or    STATE3
%   actuator that can control the heat pump.            STATE4
% The sensor senses the room from an Erlang distributed amount of time (Erlang<l=0.1 s-1,k=3>). 
% Then the CPU works for a uniform distributed amount of time (Uniform<a=10 s, b=20 s>), and after that: 
%   it returns sensing with probability p1 = 50%
%   it activates the air conditioning with probability p2 = 30%
%   turns on the heat pump with probability p3 = 20%. 
% The actuators take an exponentially distributed amount of time, respectively with rates 
%   (Exp<l=0.05 s-1>) for the air conditioning
%   (Exp<l=0.03 s-1>) for the heat pump

% • Draw a state machine based model of the system
% • Implement it in a programming language of your choice
% • Compute the probability of the system being sensing, using the CPU, actuating the air conditioning or the heat pump.
% • Determine the sensing frequency (throughput) of the system, measured in times the system enters the sensing state per time unit.

clear all;

s=1;
t=0;
i=0;
Tmax=10000;
Ts1=0;
Ts2=0;
Ts3=0;
Ts4=0;
C = 0;


while t<Tmax
    if s == 1
	   ns=2;
       dt=-(log(rand())+ log(rand())+ log(rand()))/0.1; %generating three stage Erlang dist (Erlang<l=0.1 s-1,k=3>)
       Ts1=Ts1+dt;
    end 
    if s==2
        r = rand();
        if r<0.5
               ns=1;
               C = C + 1;
        elseif r>0.5 && r<0.8
               ns=3;
        elseif r>0.8
               ns=4;
        end
        dt=10 + (20-10)*rand();  %generate Uniform distribution (Uniform<a=10 s, b=20 s>)
        Ts2=Ts2+dt;
    end
    if s==3
        ns=1;
        dt=-log(rand())/0.05;  %generate Exp distribution (Exp<l=0.05 s-1>) 
        Ts3=Ts3+dt;
        C = C + 1;
    end
    if s==4
        ns=1;
        dt=-log(rand())/0.03;  %generate Exp distribution (Exp<l=0.03 s-1>)
        Ts4=Ts4+dt;
        C = C + 1;
    end
   s=ns;
   t=t+dt;
end


% • Compute the probability of the system being sensing, using the CPU, actuating the air conditioning or the heat pump.
Ps1=Ts1/t;   %prob of being in state 1
fprintf("Probability of the system being sensing: %f\n", Ps1);
Ps2=Ts2/t;   %prob of being in state 2
fprintf("Probability of the system using the CPU: %f\n", Ps2);
Ps3=Ts3/t;   %prob of being in state 3
fprintf("Probability of the system actuating the air conditioning: %f\n", Ps3);
Ps4=Ts4/t;   %prob of being in state 4
fprintf("Probability of the system actuating the heat pump: %f\n", Ps4);


% • Determine the sensing frequency (throughput) of the system, measured in times the system enters the sensing state per time unit.
X = C / t;
fprintf("Sensing frequency (throughput) of the system: %f\n", X);
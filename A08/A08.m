%% Batch Processing with Garbage collector
% A batch processing system, requires an Exp<0.05 s-1> distributed time to prepare a new job for execution. 
% Processes run on a Java Virtual Machine, and the execution of the job requires a 
% different amount of time, depending on whether the garbage collector is running or not. 
% In particular it takes: 
% - Exp<1 s-1> when running at full speed
% - Exp<0.3 s-1> during garbage collection. 
% Garbage collection is started every Exp<0.1 s-1> and lasts Exp<0.4 s-1>.

% • Draw a state machine based model of the system
% • Implement it in a programming language of your choice
% • Compute the probability of preparing a new job, executing it at full speed, and running the task during garbage collection.
% • Determine the throughput of the system (how many new jobs per second are started)

clear all;

s=1;
t=0;
Tmax=100000;
Ts1=0;  %preparing a new job
Ts2=0;  %executing a job at full speed
Ts3=0;  %executing a job during garbage collection
j=0;    %new jobs per second


while t<Tmax
    if s==1
	   ns=2;
       dt= -log(rand())/0.05; %Exp<0.05 s-1> distributed time to prepare a new job for execution
       Ts1=Ts1+dt;
       j=j+1;
    end 
    if s==2
        dtFullSpeed=-log(rand())/1;   %job execution takes Exp<1 s-1> when running at full speed
        dtGCstart=-log(rand())/0.1;   %Garbage collection is started every Exp<0.1 s-1>
        dt=min(dtFullSpeed,dtGCstart);
        Ts2=Ts2+dt;
        if dtFullSpeed<dtGCstart
          ns=1;
        else 
          ns=3;
        end
    end
     if s==3
        dtSlowSpeed=-log(rand())/0.3;   %job execution takes Exp<0.3 s-1> during garbage collection
        dtGCend=-log(rand())/0.4;   %Garbage collection lasts Exp<0.4 s-1>
        dt=min(dtSlowSpeed,dtGCend);
        Ts3=Ts3+dt;
        if dtSlowSpeed<dtGCend
          ns=1;
        else 
          ns=2;
        end
     end
   s=ns;
   t=t+dt;
end

Ps1=Ts1/t;   %probability of preparing a new job
fprintf("Probability of preparing a new job: %f\n", Ps1);
Ps2=Ts2/t;   %probability of executing a job at full speed
fprintf("Probability of executing a job at full speed: %f\n", Ps2);
Ps3=Ts3/t;   %probability of executing a job during garbage collection
fprintf("Probability of executing a job during garbage collection: %f\n", Ps3);

X = j/Tmax;
fprintf("Throughput of the system (how many new jobs per second are started): %f\n", X);
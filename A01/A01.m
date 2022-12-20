clear all;

filename1 = 'apache1.log';
filename2 = 'apache2.log';

data = readtable(filename1, 'FileType', 'text'); % Replace "filename1" with "filename2" to compute Basic Performance Indices and Workloads of apache2.log
A_data = data(:,4);
A_str = table2array(A_data);

A_split = split(A_str(:,end),'[');
A_str = A_split(:,2);
A = datetime(A_str,'InputFormat','dd/MMM/uuuu:HH:mm:ss.SSS', 'Format', 'dd/MMM/uuuu:HH:mm:ss.SSS');

S_data = data(:, end);
S = table2array(S_data);

nA = size(data, 1);

At = milliseconds(A(:) - A(1)); % Delta t of arrivals wrt first arrival set to 0

ASC = [At, S, zeros(nA,1)];

ASC(1,3) = sum(ASC(1,1:2));

for i = 1:nA-1
    if ASC(i+1, 1) <= ASC(i,3)
       ASC(i+1, 3) = ASC(i,3) + ASC(i+1, 2);   
    else
        ASC(i+1, 3) = ASC(i+1,1) + ASC(i+1, 2);
    end

end



T_ms = ASC(nA,3); % T in ms
T = T_ms/1000; % T in s

lambda = nA/T; % Arrival rate

nC = nA;
X = nC/T; % Throughput

interArr = 1/lambda;   % Average inter-arrival time

B = sum(ASC(:,2))/1000; % Busy time

U = B / T; % Utilization rate in s

Rt = (ASC(:,3) - ASC(:,1))/1000; % Rt in s

W = sum(Rt); % Area of the difference between arrivals and departures

S = B / nC; % Average service time

N = W / T;  % Average number of jobs

R = N / X;  % Average response time


evs = [ASC(:,1), ones(nA, 1), zeros(nA, 6); 
       ASC(:,3), -ones(nC,1), zeros(nC, 6)];
evs = sortrows(evs, 1);
evs(:,3) = cumsum(evs(:,2)); % NofT
evs(1:end-1, 4) = evs(2:end,1) - evs(1:end-1,1); % dT
evs(:,5) = (evs(:,3) == 0) .* evs(:,4);
evs(:,6) = (evs(:,3) == 1) .* evs(:,4);
evs(:,7) = (evs(:,3) == 2) .* evs(:,4);
evs(:,8) = (evs(:,3) == 3) .* evs(:,4);



Y1 = sum(evs(:,5)); % Fraction of time the system has been with 0 jobs
p1 = Y1 / T_ms;

Y2 = sum(evs(:,6));
p2 = Y2 / T_ms;

Y3 = sum(evs(:,7));
p3 = Y3 / T_ms;

Y4 = sum(evs(:,8));
p4 = Y4 / T_ms;

p5 = sum(Rt < 1) / nC;

p6 = sum(Rt < 5) / nC;

p7 = sum(Rt < 10) / nC;

fprintf(1, "Basic Performance Indices and Workloads \n\n");

fprintf(1, "Arrival rate: %g\n", lambda);
fprintf(1, "Throughput: %g\n", X);
fprintf(1, "Average inter arrival time: %g\n", interArr);
fprintf(1, "Busy time: %g\n", B);
fprintf(1, "Utilization: %g\n", U);
fprintf(1, "Residence time (area of the difference between arrivals and departures): %g\n", W);
fprintf(1, "Average service time: %g\n", S);
fprintf(1, "Average number of jobs: %g\n", N);
fprintf(1, "Average response time: %g\n", R);

fprintf(1, "Probability of having m jobs in the web server with m = 0: %g\n", p1);
fprintf(1, "Probability of having m jobs in the web server with m = 1: %g\n", p2);
fprintf(1, "Probability of having m jobs in the web server with m = 2: %g\n", p3);
fprintf(1, "Probability of having m jobs in the web server with m = 3: %g\n", p4);
fprintf(1, "Probability of having a response time less than tao = 1 s: %g\n", p5);
fprintf(1, "Probability of having a response time less than tao = 5 s: %g\n", p6);
fprintf(1, "Probability of having a response time less than tao = 10 s: %g\n", p7);

%% Performance indices of the Intranet of a school
% The intranet of a school, with N = 530 students, is composed by three servers: 
%   - the Moodle server
%   - the file server
%   - the DB server. 
% Their average service times are respectively: 
%   - Sm = 80ms
%   - Sf = 120ms 
%   - Sd = 11ms. 
% The three stations are also characterized by the following visits: 
%   - vm = 1
%   - vf = 0.75
%   - vd = 10. 
% Knowing that the think time of the students is Z = 2 min, compute:
%   1. The demand of the three station
%   2. The system throughput
%   3. The utilization of the three stations
%   4. The system response time
%   5. The average number of students "not thinking"

% number of students 
N = 530;
% average Service times [ms --> s]
S = [80, 120, 11];
S = S * 10^-3;
% visists
v = [1, 0.75, 10];
% think time [min --> s]
Z = 2 * 60;

%   1. The demand of the three station
D = S.*v;
disp("D")
disp(D)

% Mean Value Analysis
Qk = zeros(1,3);
Rk= zeros(1,3);

for i = 1:N
    for j = 1:3
       Rk(j) = D(j)*(1+Qk(1,j)); 
    end
    X = i/(Z+sum(Rk)); 
    for k = 1:3
        Qk(k) = X*Rk(k); 
    end
end

%   2. The system throughput
disp("X")
disp(X)

%   3. The utilization of the three stations
Xk = v * X;
Uk = Xk .* S;
disp("Uk:")
disp(Uk)

%   4. The system response time
R = sum(Rk);
disp("R")
disp(R)

%   5. The average number of students "not thinking"
N_notThinking = X * R;
disp("N_notThinking")
disp(N_notThinking)








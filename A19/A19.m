%% Response time of a storage system
% A small company mainly uses:
%   - a Costumer Relation Management (CRM – resource 1) software, 
%   - and a File Sharing (FS – resource 2) service. 
% Each of them is hosted on a separate server, which works in processor sharing, and serves two type of users: 
%   - employees 
%   - and providers. 
% Requests of both type of users arrive according to Poisson processes respectively of rates:
%   - lambdaE = 5 req./s. 
%   - lambdaP = 10 req./s.
% Since the CRM uses resources on the FS, requests need always to be handled by both servers. 
% The demand at each server for the two types of users have been measured as follows:
%   - D1E = 50 ms.      D1P = 60 ms.
%   - D2E = 100 ms.     D2P = 40 ms.
% Determine (using a programming language, and NOT a tool like a JMT):
%   1. The utilization of the two servers
%   2. The residence time of the two servers
%   3. The system response time

lambdaE = 5;   %req./s
lambdaP = 10;  %req./s
D1E = 50/1000;  %ms --> s      
D1P = 60/1000;  %ms --> s
D2E = 100/1000; %ms --> s   
D2P = 40/1000;  %ms --> s

%   1. The utilization of the two servers
U1 = lambdaE * D1E + lambdaP * D1P;
U2 = lambdaE * D2E + lambdaP * D2P;
disp("U1")
disp(U1)
disp("U2")
disp(U2)
    
%   2. The residence time of the two servers
R1E = D1E / (1 - U1);
R1P = D1P / (1 - U1);
R2E = D2E / (1 - U2);
R2P = D2P / (1 - U2);
R1 = lambdaE / (lambdaE + lambdaP) * R1E + lambdaP / (lambdaE + lambdaP) * R1P;
R2 = lambdaE / (lambdaE + lambdaP) * R2E + lambdaP / (lambdaE + lambdaP) * R2P;
disp("R1")
disp(R1)
disp("R2")
disp(R2)

%   3. The system response time
R = R1 + R2;
disp("R")
disp(R)
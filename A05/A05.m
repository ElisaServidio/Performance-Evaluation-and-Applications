clear all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% For each distribution, generate N = 500 samples using the techniques seen during the course, taking samples from the uniform distribution from the cells of the enclosed file. 
% In particular, use the first column for selecting branches for discrete distributions (cases 2, 4 and 6), and the other two columns for generating samples of continuous distributions.
% Hint: one single random number is required for cases 1, 2 and 3; two values are required for cases 4 and 5 – respectively columns one and two, and columns two and three; 
% case 6 requires the use of all three values – one for selecting the branch, the one in column two for the first branch, and the ones in columns two and three for the second branch.

% Plot a figure where you compare the CDF from the generated samples, with the analytical one using the given parameters.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

filename = 'random.csv';

%import data
data = readtable(filename); 
X = table2array(data);
N = 500;


% 1. A continuous uniform distribution between [5, 15]
a = 5;
b = 15;
res1 = a + (b-a) * (X(:,2));

xunif = [5:15];              % values and probabilities to plot analytical CDF
punif = (xunif-a)/(b-a);

% 2. A discrete distribution with Value Probability
%                                   5       0.3
%                                  10       0.4
%                                  15       0.3
x = [5, 10, 15];
p2 = [0.3, 0.4, 0.3];
C2 = cumsum(p2);

for k = 1:N
    r = X(k,1);
    for i = 1:3
        if r <= C2(1,i)
            res2(k,1) = x(1,i);
            break
        end
    end
end

xdiscr = [0:15];              % values and probabilities to plot analytical CDF
pdiscr = zeros(size(xdiscr)); %descrete CDF = P(X <= x) (stairs plot)
for i = 1:16
    if i-1 >= 5 && i-1 < 10
        pdiscr(1,i) = C2(1,1);
    end
    if i-1 >= 10 && i-1 < 15
        pdiscr(1,i) = C2(1,2);
    end
    if i-1 >= 15
        pdiscr(1,i) = C2(1,3);
    end
end

% 3. An exponential distribution with average 10 
lambda = 1/10; 
res3 = - log(X(:,2)) / lambda;

xexp = [0:80];                  % values and probabilities to plot analytical CDF
pexp = 1-exp(-lambda * xexp);

% 4. An hyper-exponential distribution with two stages, characterized by (l1 = 0.05, l2 = 0.175, p1 = 0.3)
lhyper = [0.05, 0.175];
p1hyper = [0.3, 0.7];
C3 = cumsum(p1hyper);

for k = 1:N
    r = X(k,1);
    for i = 1:2
        if r <= C3(1,i)
            res4(k,1) = - log(X(k,2)) / lhyper(1,i);
            break
        end
    end
end

xhypexp = [0:80];                   % values and probabilities to plot analytical CDF
phypexp = 1 - (p1hyper(1,1) * exp(-lhyper(1,1) * xhypexp) + (1 - p1hyper(1,1)) * exp(-lhyper(1,2) * xhypexp));

% 5. An hypo-exponential distribution with two stages characterized by (l1 = 0.25, l2 = 0.16667)
l1hypo = 0.25;
l2hypo = 0.16667;
res5 = - log(X(:,2)) / l1hypo - log(X(:,3)) / l2hypo;

xhypoexp = [0:60];              % values and probabilities to plot analytical CDF
phypoexp = 1 - (l2hypo * exp(-l1hypo * xhypoexp)) / (l2hypo - l1hypo) + (l1hypo * exp(-l2hypo * xhypoexp)) / (l2hypo - l1hypo);

% 6. An Hyper-Erlang characterized by the following branches, Num. stages (k) Rate (l) Probability
%                                                                    1         0.05        0.3
%                                                                    2         0.35        0.7
K = [1, 2];
l = [0.05, 0.35];
p6 = [0.3, 0.7];
C6 = cumsum(p6);

r = X(:,1);
res6 = (r < p6(1,1)) .* (-log(X(:,2))/l(1,1)) + (r >= p6(1,1)) .* (-log(X(:,2).*X(:,3))/l(1,2));

xhyperErlang = [0:80];                 % values and probabilities to plot analytical CDF
phyperErlang = 1- p6(1,1) .* exp(-l(1,1)*xhyperErlang) - p6(1,2) .* (exp(-l(1,2)*xhyperErlang) + l(1,2)*xhyperErlang.*exp(-l(1,2)*xhyperErlang));

% Plot the distributions 
tl = tiledlayout(2,3);
nexttile
plot(sort(res1), [1:N]/N, "-", xunif, punif, "-", "LineWidth",1)
title('Continuous uniform distribution')
legend('CDF from generated samples','Analytical CDF','Location','southeast')


%%% Compute probability for plotting generated descrete CDF %%%
% In alternative, it could have been done automatically by using 
% cdfplot() function which automatically computes CDF of a dataset
cont5 = 0;
cont10 = 0;
cont15 = 0;
for i = 1:N
    if res2(i,1) == 5
        cont5 = cont5 +1;
    end
    if res2(i,1) == 10
        cont10 = cont10 +1;
    end
    if res2(i,1) == 15
        cont15 = cont15 +1;
    end
end

prob2 = [cont5/N,cont10/N,cont15/N];
cum_prob2 = cumsum(prob2);
gen_pdiscr = zeros(size(xdiscr)); %descrete CDF = P(X <= x) (stairs plot)
for i = 1:16
    if i-1 >= 5 && i-1 < 10
        gen_pdiscr(1,i) = cum_prob2(1,1);
    end
    if i-1 >= 10 && i-1 < 15
        gen_pdiscr(1,i) = cum_prob2(1,2);
    end
    if i-1 >= 15
        gen_pdiscr(1,i) = cum_prob2(1,3);
    end
end

nexttile
stairs(xdiscr, gen_pdiscr,"LineWidth",1) 
hold on
stairs(xdiscr, pdiscr,"LineWidth",1)
axis([5 15 0 1])
hold off
title('Discrete distribution')
legend('CDF from generated samples','Analytical CDF','Location','southeast')


nexttile
plot(sort(res3), [1:N]/N, "-", xexp, pexp, "-","LineWidth",1)
title('Exponential distribution')
legend('CDF from generated samples','Analytical CDF','Location','southeast')


nexttile
plot(sort(res4), [1:N]/N, "-", xhypexp, phypexp, "-","LineWidth",1)
title('Hyper-exponential distribution')
legend('CDF from generated samples','Analytical CDF','Location','southeast')


nexttile
plot(sort(res5), [1:N]/N, "-", xhypoexp, phypoexp, "-","LineWidth",1)
title('Hypo-exponential distribution')
legend('CDF from generated samples','Analytical CDF','Location','southeast')


nexttile
plot(sort(res6), [1:N]/N, "-", xhyperErlang, phyperErlang, "-","LineWidth",1)
title('Hyper-Erlang distribution')
legend('CDF from generated samples','Analytical CDF','Location','southeast')


%Print the first 10 results of each distribution
fprintf(1, "Res1: \n")
disp(res1(1:10));
fprintf("\n")

fprintf(1, "Res2: \n")
disp(res2(1:10));
fprintf("\n")

fprintf(1, "Res3: \n")
disp(res3(1:10));
fprintf("\n")

fprintf(1, "Res4: \n")
disp(res4(1:10));
fprintf("\n")

fprintf(1, "Res5: \n")
disp(res5(1:10));
fprintf("\n")

fprintf(1, "Res6: \n")
disp(res6(1:10));
fprintf("\n")

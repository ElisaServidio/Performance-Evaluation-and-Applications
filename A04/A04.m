
clear all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% - Fit the data against the uniform, the exponential, the two stages hyper-exponential and the two stage hypo-exponential distributions using the method of moments

% - Fit the data against the exponential, the two stages hyper-exponential and the two stage hypo-exponential distributions using the method of maximum likelihood

% - With a plot, graphically compare the CDF of the considered distributions, using the parameters obtained by both fitting procedures. 
%   Make two plots per trace, the first showing the results from the method of moments, and the second for the maximum likelihood. 
%   Each plot, must show the results of the corresponding fitting, and include the approximated CDF obtained directly from the dataset.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

filename = 'Traces.csv';

%import data
data = readtable(filename); 
X = table2array(data);
Y = sort(X);

% METHOD OF MOMENTS

M1 = sum(X)  / 1000;
M2 = sum(X.^2) / 1000;
M3 = sum(X.^3) / 1000;

cv = sqrt(M2 - M1.^2) ./ M1;

% Exponential distribution

lambda_exp = 1 ./ M1;

t = [0:0.1:45]; 

y1exp = 1- exp(-lambda_exp(1,1) * t);
y2exp = 1- exp(-lambda_exp(1,2) * t);
y3exp = 1- exp(-lambda_exp(1,3) * t);

% Uniform distribution

a_unif = M1 - 0.5 * sqrt(12 * (M2 - M1.^2));
b_unif = M1 + 0.5 * sqrt(12 * (M2 - M1.^2));

y1unif = max(0, min(1, (t - a_unif(1,1)) / (b_unif(1,1) - a_unif(1,1))));
y2unif = max(0, min(1, (t - a_unif(1,2)) / (b_unif(1,2) - a_unif(1,2))));
y3unif = max(0, min(1, (t - a_unif(1,3)) / (b_unif(1,3) - a_unif(1,3))));

% Two stages hyper-exponential distribution
% An Hyper-exponential will never have a c.v. < 1

pars_hyper = [0.8/M1(1,1), 1.2/M1(1,1), 0.4]; % it is important to give meaningful starting points for the parameters

%Create a procedure that returns the moments of the distribution
HypExpMoments = @(p) [p(1,3)/p(1,1) + (1-p(1,3))/p(1,2), 2*(p(1,3)/p(1,1)^2 + (1-p(1,3))/p(1,2)^2), 6*(p(1,3)/p(1,1)^3 + (1-p(1,3))/p(1,2)^3)];
HypExpMoments(pars_hyper);

F1 = @(p) (HypExpMoments(p) ./ [M1(1,1), M2(1,1), M3(1,1)] - 1);
F1(pars_hyper);
F2 = @(p) (HypExpMoments(p) ./ [M1(1,2), M2(1,2), M3(1,2)] - 1);
F2(pars_hyper);
F3 = @(p) (HypExpMoments(p) ./ [M1(1,3), M2(1,3), M3(1,3)] - 1);
F3(pars_hyper);

res1HypExp = fsolve(F1, pars_hyper); 
res2HypExp = fsolve(F2, pars_hyper); 
res3HypExp = fsolve(F3, pars_hyper); 
 
HypExpCDF = @(p, t) 1 - p(1,3)*exp(-p(1,1)*t) - p(1,3)*exp(-p(1,2)*t);

y1hypexp = HypExpCDF(res1HypExp, t);
y2hypexp = HypExpCDF(res2HypExp, t);
y3hypexp = HypExpCDF(res3HypExp, t);

% Two stage hypo-exponential distribution
% An Hypo-exponential cannot have c.v. > 1

pars_hypo = [1/(0.3*M1(1,1)), 1/(0.7*M1(1,1))]; % it is important to give meaningful starting points for the parameters

%Create a procedure that returns the moments of the distribution
HypoExpMoments = @(p) [1/(p(1,1) - p(1,2))*(p(1,1)/p(1,2)-p(1,2)/p(1,1)),2*(1/p(1,1)^2+1/(p(1,1)*p(1,2))+1/p(1,2)^2)];
HypoExpMoments(pars_hypo);

F1 = @(p) (HypoExpMoments(p) ./ [M1(1,1), M2(1,1)] - 1);
F1(pars_hypo);
F2 = @(p) (HypoExpMoments(p) ./ [M1(1,2), M2(1,2)] - 1);
F2(pars_hypo);
F3 = @(p) (HypoExpMoments(p) ./ [M1(1,3), M2(1,3)] - 1);
F3(pars_hypo);

res1HypoExp = fsolve(F1, pars_hypo); 
res2HypoExp = fsolve(F2, pars_hypo); 
res3HypoExp = fsolve(F3, pars_hypo);  

HypoExpCDF = @(p, t) 1 - p(1,2)*exp(-p(1,1)*t)/(p(1,2)-p(1,1)) + p(1,1)*exp(-p(1,2)*t)/(p(1,2)-p(1,1));

y1hypoexp = HypoExpCDF(res1HypoExp, t);
y2hypoexp = HypoExpCDF(res2HypoExp, t);
y3hypoexp = HypoExpCDF(res3HypoExp, t);

% METHOD OF MAXIMUM LIKELIHOOD

% Two stages hyper-exponential distribution

HypeExpPDF = @(x, p1, l1, l2) p1*l1*exp(-x*l1)+(1-p1)*l2*exp(-x*l2);
res_1HypExp = mle(X(:,1), 'pdf', HypeExpPDF, 'Start', [0.4, 0.8/M1(1,1), 1.2/M1(1,1)], 'LowerBound', [0.0001,0.0001, 0.0001], 'UpperBound', [1,Inf,Inf]);
res_2HypExp = mle(X(:,2), 'pdf', HypeExpPDF, 'Start', [0.4, 0.8/M1(1,2), 1.2/M1(1,2)], 'LowerBound', [0.0001,0.0001, 0.0001], 'UpperBound', [1,Inf,Inf]);
res_3HypExp = mle(X(:,3), 'pdf', HypeExpPDF, 'Start', [0.4, 0.8/M1(1,3), 1.2/M1(1,3)], 'LowerBound', [0.0001,0.0001, 0.0001], 'UpperBound', [1,Inf,Inf]);

y1hypexp_mle = HypExpCDF(res_1HypExp, t);
y2hypexp_mle = HypExpCDF(res_2HypExp, t);
y3hypexp_mle = HypExpCDF(res_3HypExp, t);

% Two stage hypo-exponential distribution

HypoExpPDF = @(x, l1, l2) (l1*l2 / (l1-l2))*(exp(-x*l2)-exp(-x*l1));
res_1HypoExp = mle(X(:,1), 'pdf', HypoExpPDF, 'Start', [1/(0.3*M1(1,1)), 1/(0.7*M1(1,1))], 'LowerBound', [0.0001,0.0001], 'UpperBound', [Inf,Inf]);
res_2HypoExp = mle(X(:,2), 'pdf', HypoExpPDF, 'Start', [1/(0.3*M1(1,2)), 1/(0.7*M1(1,2))], 'LowerBound', [0.0001,0.0001], 'UpperBound', [Inf,Inf]);
res_3HypoExp = mle(X(:,3), 'pdf', HypoExpPDF, 'Start', [1/(0.3*M1(1,3)), 1/(0.7*M1(1,3))], 'LowerBound', [0.0001,0.0001], 'UpperBound', [Inf,Inf]);

y1hypoexp_mle = HypoExpCDF(res_1HypoExp, t);
y2hypoexp_mle = HypoExpCDF(res_2HypoExp, t);
y3hypoexp_mle = HypoExpCDF(res_3HypoExp, t);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp("Column A");

disp("METHOD OF MOMENTS:");
fprintf(1, "1st moment: %g\n", M1(1,1));
fprintf(1, "2nd moment: %g\n", M2(1,1));
fprintf(1, "3rd moment: %g\n", M3(1,1));
fprintf(1, "Coefficient of variation: %g\n", cv(1,1));
fprintf(1, "Uniform Left bound a: %g\n", a_unif(1,1));
fprintf(1, "Uniform Right bound b: %g\n", b_unif(1,1));
fprintf(1, "Exponential rate lambda: %g\n", lambda_exp(1,1));
if (cv(1,1) < 1)
    disp ("No real solution for Hyper-exponential because cv < 1");
else
    fprintf(1, "Hyper Exponential lambda1: %g\n", res1HypExp(1,1));
    fprintf(1, "Hyper Exponential lambda2: %g\n", res1HypExp(1,2));
    fprintf(1, "Hyper Exponential p1: %g\n", res1HypExp(1,3));
end
if (cv(1,1) > 1)
    disp ("No real solution for Hypo-exponential because cv > 1");
    fprintf(1, '\n');
else
    fprintf(1, "Hypo Exponential lambda1: %g\n", res1HypoExp(1,1));
    fprintf(1, "Hypo Exponential lambda2: %g\n\n", res1HypoExp(1,2));
end
 
disp("METHOD OF MAXIMUM LIKELIHOOD:");
fprintf(1, "Exponential rate lambda: %g\n", lambda_exp(1,1));
if (cv(1,1) < 1)
    disp ("No real solution for Hyper-exponential because cv < 1");
else
    fprintf(1, "Hyper Exponential lambda1: %g\n", res_1HypExp(1,2));
    fprintf(1, "Hyper Exponential lambda2: %g\n", res_1HypExp(1,3));
    fprintf(1, "Hyper Exponential p1: %g\n", res_1HypExp(1,1));
end

if (cv(1,1) > 1)
    disp ("No real solution for Hypo-exponential because cv > 1");
    fprintf(1, '\n');
else
    fprintf(1, "Hypo Exponential lambda1: %g\n", res_1HypoExp(1,1));
    fprintf(1, "Hypo Exponential lambda2: %g\n\n", res_1HypoExp(1,2));
end

disp("Column B");

disp("METHOD OF MOMENTS:");
fprintf(1, "1st moment: %g\n", M1(1,2));
fprintf(1, "2nd moment: %g\n", M2(1,2));
fprintf(1, "3rd moment: %g\n", M3(1,2));
fprintf(1, "Coefficient of variation: %g\n", cv(1,2));
fprintf(1, "Uniform Left bound a: %g\n", a_unif(1,2));
fprintf(1, "Uniform Right bound b: %g\n", b_unif(1,2));
fprintf(1, "Exponential rate lambda: %g\n", lambda_exp(1,2));
if (cv(1,2) < 1)
    disp ("No real solution for Hyper-exponential because cv < 1");
else
    fprintf(1, "Hyper Exponential lambda1: %g\n", res2HypExp(1,1));
    fprintf(1, "Hyper Exponential lambda2: %g\n", res2HypExp(1,2));
    fprintf(1, "Hyper Exponential p1: %g\n", res2HypExp(1,3));
end
if (cv(1,2) > 1)
    disp ("No real solution for Hypo-exponential because cv > 1");
    fprintf(1, '\n');
else
    fprintf(1, "Hypo Exponential lambda1: %g\n", res2HypoExp(1,1));
    fprintf(1, "Hypo Exponential lambda2: %g\n\n", res2HypoExp(1,2));
end
 
disp("METHOD OF MAXIMUM LIKELIHOOD:");
fprintf(1, "Exponential rate lambda: %g\n", lambda_exp(1,2));
if (cv(1,2) < 1)
    disp ("No real solution for Hyper-exponential because cv < 1");
else
    fprintf(1, "Hyper Exponential lambda1: %g\n", res_2HypExp(1,2));
    fprintf(1, "Hyper Exponential lambda2: %g\n", res_2HypExp(1,3));
    fprintf(1, "Hyper Exponential p1: %g\n", res_2HypExp(1,1));
end

if (cv(1,2) > 1)
    disp ("No real solution for Hypo-exponential because cv > 1");
    fprintf(1, '\n');
else
    fprintf(1, "Hypo Exponential lambda1: %g\n", res_2HypoExp(1,1));
    fprintf(1, "Hypo Exponential lambda2: %g\n\n", res_2HypoExp(1,2));
end

disp("Column C");

disp("METHOD OF MOMENTS:");
fprintf(1, "1st moment: %g\n", M1(1,3));
fprintf(1, "2nd moment: %g\n", M2(1,3));
fprintf(1, "3rd moment: %g\n", M3(1,3));
fprintf(1, "Coefficient of variation: %g\n", cv(1,3));
fprintf(1, "Uniform Left bound a: %g\n", a_unif(1,3));
fprintf(1, "Uniform Right bound b: %g\n", b_unif(1,3));
fprintf(1, "Exponential rate lambda: %g\n", lambda_exp(1,3));
if (cv(1,3) < 1)
    disp ("No real solution for Hyper-exponential because cv < 1");
else
    fprintf(1, "Hyper Exponential lambda1: %g\n", res3HypExp(1,1));
    fprintf(1, "Hyper Exponential lambda2: %g\n", res3HypExp(1,2));
    fprintf(1, "Hyper Exponential p1: %g\n", res3HypExp(1,3));
end
if (cv(1,3) > 1)
    disp ("No real solution for Hypo-exponential because cv > 1");
    fprintf(1, '\n');
else
    fprintf(1, "Hypo Exponential lambda1: %g\n", res3HypoExp(1,1));
    fprintf(1, "Hypo Exponential lambda2: %g\n\n", res3HypoExp(1,2));
end
 
disp("METHOD OF MAXIMUM LIKELIHOOD:");
fprintf(1, "Exponential rate lambda: %g\n", lambda_exp(1,3));
if (cv(1,3) < 1)
    disp ("No real solution for Hyper-exponential because cv < 1");
else
    fprintf(1, "Hyper Exponential lambda1: %g\n", res_3HypExp(1,2));
    fprintf(1, "Hyper Exponential lambda2: %g\n", res_3HypExp(1,3));
    fprintf(1, "Hyper Exponential p1: %g\n", res_3HypExp(1,1));
end

if (cv(1,3) > 1)
    disp ("No real solution for Hypo-exponential because cv > 1");
    fprintf(1, '\n');
else
    fprintf(1, "Hypo Exponential lambda1: %g\n", res_3HypoExp(1,1));
    fprintf(1, "Hypo Exponential lambda2: %g\n\n", res_3HypoExp(1,2));
end


% Figure for the Method of moments
% It must show uniform, exponential, the two stages hyper-exponential and the 
% two stage hypo-exponential distributions using the method of moments
% and include the approximated CDF obtained directly from the dataset.
tl = tiledlayout(2,3);
nexttile
plot(Y(:,1), [1:1000]/1000, "+", 'MarkerSize', 6, "Color", "cyan")
hold on
plot(t, y1exp, "-", 'LineWidth', 3, "Color", "magenta")
hold on
plot(t, y1unif, "-", 'LineWidth', 1,"Color", "blue")
hold on
plot(t, y1hypoexp, "-", 'LineWidth', 0.7, "Color", "black")
hold off
title('Distribution Column A with Method of moments')
legend('dataset', 'exponential', 'uniform', 'hypo-exponential')
axis([0 45 0 1])
nexttile
plot(Y(:,2), [1:1000]/1000, "+", 'MarkerSize', 6, "Color", "cyan")
hold on
plot(t, y2exp, "-", 'LineWidth', 1, "Color", "magenta")
hold on
plot(t, y2unif, "-", 'LineWidth', 1,"Color", "blue")
hold on
plot(t, y2hypexp, "-", 'LineWidth', 1, "Color", "black")
hold off
title('Distribution Column B with Method of moments')
legend('dataset', 'exponential', 'uniform', 'hyper-exponential')
axis([0 45 0 1])
nexttile
plot(Y(:,3), [1:1000]/1000, "+", 'MarkerSize', 6, "Color", "cyan")
hold on
plot(t, y3exp, "-", 'LineWidth', 1, "Color", "magenta")
hold on
plot(t, y3unif, "-", 'LineWidth', 1,"Color", "blue")
hold on
plot(t, y3hypoexp, "-", 'LineWidth', 1, "Color", "black")
hold off
title('Distribution Column C with Method of moments')
legend('dataset', 'exponential', 'uniform', 'hypo-exponential')
axis([0 45 0 1])


% Figure for the Maximum likelihood
% It must show the exponential, the two stages hyper-exponential and the 
% two stage hypo-exponential distributions using the method of maximum likelihood
% and include the approximated CDF obtained directly from the dataset.
nexttile
plot(Y(:,1), [1:1000]/1000, "+", 'MarkerSize', 6, "Color", "cyan")
hold on
plot(t, y1exp, "-", 'LineWidth', 3, "Color", "magenta")
hold on
plot(t, y1hypoexp_mle, "-", 'LineWidth', 0.7, "Color", "black")
hold off
title('Distribution Column A with Maximum likelihood')
legend('dataset', 'exponential', 'hypo-exponential')
axis([0 45 0 1])
nexttile
plot(Y(:,2), [1:1000]/1000, "+", 'MarkerSize', 6, "Color", "cyan")
hold on
plot(t, y2exp, "-", 'LineWidth', 1, "Color", "magenta")
hold on
plot(t, y2hypexp_mle, "-", 'LineWidth', 1, "Color", "black")
hold off
title('Distribution Column B with Maximum likelihood')
legend('dataset', 'exponential', 'hyper-exponential')
axis([0 45 0 1])
nexttile
plot(Y(:,3), [1:1000]/1000, "+", 'MarkerSize', 6, "Color", "cyan")
hold on
plot(t, y3exp, "-", 'LineWidth', 1, "Color", "magenta")
hold on
plot(t, y3hypoexp_mle, "-", 'LineWidth', 1, "Color", "black")
hold off
title('Distribution Column C with Maximum likelihood')
legend('dataset', 'exponential', 'hypo-exponential')
axis([0 45 0 1])
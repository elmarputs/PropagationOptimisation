%% Generate composite design

% Number of control variables
N = 6;
% Create variable settings matrix
variableSettings = ccdesign(N);

% Export matrix to file
writematrix(variableSettings);

%% Import central composite design results

ccdResults = readmatrix('ccdFitness.txt');

decisionVariables = ccdResults(:, 2:7);
fitness = ccdResults(:, 8:9);
fitness(:,1) = fitness(:,1)./1e9;
fitness(:,2) = fitness(:,2)./1e6;

rstool(decisionVariables, fitness, 'quadratic')

%% Compute R² (goodness of fit)

Sr1 = sumsqr(residuals(:,1));
Sr2 = sumsqr(residuals(:,2));

dev1 = fitness(:,1) - mean(fitness(:,1));
dev2 = fitness(:,2) - mean(fitness(:,2));
St1 = sumsqr(dev1);
St2 = sumsqr(dev2);

R1 = 1 - Sr1/St1;
R2 = 1 - Sr2/St2;

N = length(fitness(:,1));
p = length(beta);
Ra1 = 1-(1-R1)*((N-1)/(N-p));
Ra2 = 1-(1-R2)*((N-1)/(N-p));

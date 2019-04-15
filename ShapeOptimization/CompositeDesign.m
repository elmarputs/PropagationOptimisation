%% Generate composite design

% Number of control variables
N = 6;
% Create variable settings matrix
variableSettings = ccdesign(N);

% Export matrix to file
writematrix(variableSettings);

%% Import central composite design results

ccdResults = readmatrix('ccdFitness.txt');

xValues = ccdResults(:, 2);
yValues = ccdResults(:, 3);
x = linspace(min(xValues), max(xValues), 50);
y = linspace(min(yValues), max(yValues), 50);
fitness = ccdResults(:, 4);
[X, Y] = meshgrid(x, y);
Z = griddata(xValues, yValues, fitness, X, Y);

figure;
surfc(X, Y, Z);

%% Construct response surface

% Perform least squares regression
A = [ones(length(xValues),1) xValues yValues xValues.^2 yValues.^2 xValues.*yValues.^2 xValues.^2.*yValues xValues.^3.*yValues xValues.*yValues.^3 xValues.^4 yValues.^4];
b = lsqr(A, fitness);

i = 1;
for xV = x
    j = 1;
    for yV = y
        f(i,j) = [1 xV yV xV^2 yV^2 xV*yV^2 xV^2*yV xV^3*yV xV*yV^3 xV^4 yV^4] * b;
        j = j + 1;
    end
    i = i + 1;
end

figure;
surfc(X, Y, f);

% Compute goodness of fit
fitnessMean = mean(fitness);
fitValues = A*b;
Sr = sumsqr(fitValues - fitnessMean);
St = sumsqr(fitness - fitnessMean);
R_sqr = Sr/St;
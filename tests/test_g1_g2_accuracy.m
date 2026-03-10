%% TEST: g1 and g2 Function Accuracy
%
%  Validates the closed-form approximation functions g1() and g2() against
%  Monte Carlo simulation of E[log(X)] and E[log(X+b)] for a non-central
%  chi-squared RV X ~ ncx2(2, lambda).
%
%  Expected: approximation accuracy < 1% for R >= 10 terms.
%
% -------------------------------------------------------------------------

clc; clear;

addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'src', 'analytics'));

fprintf('=== Testing g1() accuracy ===\n');
lambdaVec = [0, 1, 3, 5, 10, 20];
N = 1e6;
R = 15;

fprintf('%-10s %-15s %-15s %-10s\n', 'lambda', 'g1 Approx', 'MC Exact', 'Rel Err %');
fprintf('%s\n', repmat('-', 1, 55));
for i = 1:length(lambdaVec)
    lam = lambdaVec(i);
    X   = ncx2rnd(2, lam, [1, N]);
    mc  = mean(log(X));
    approx = g1(lam, R);
    err = abs(approx - mc) / abs(mc) * 100;
    fprintf('%-10.1f %-15.6f %-15.6f %-10.3f\n', lam, approx, mc, err);
end

fprintf('\n=== Testing g2() accuracy ===\n');
b_vals = [0.1, 0.5, 1.0, 2.0];
lam    = 5;
fprintf('%-10s %-15s %-15s %-10s\n', 'b', 'g2 Approx', 'MC Exact', 'Rel Err %');
fprintf('%s\n', repmat('-', 1, 55));
for i = 1:length(b_vals)
    b   = b_vals(i);
    X   = ncx2rnd(2, lam, [1, N]);
    mc  = mean(log(X + b));
    approx = g2(lam, b, R);
    err = abs(approx - mc) / abs(mc) * 100;
    fprintf('%-10.3f %-15.6f %-15.6f %-10.3f\n', b, approx, mc, err);
end

fprintf('\nAll tests completed.\n');

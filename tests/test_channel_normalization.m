%% TEST: Channel Normalization
%
%  Verifies that the Rician channel generation correctly achieves:
%    E[X] = L_au  (large-scale path loss)
%  after normalization.
%
%  Also verifies Rayleigh channel mean.
%
% -------------------------------------------------------------------------

clc; clear;

fprintf('=== Testing Channel Normalization ===\n\n');

% System parameters (simple test case)
D = 10; H = 1.5;
t_A = [0, 0, 0]; t_U = [D/5, 0, H];
d_au = norm(t_A - t_U);
alpha_L = 2; alpha_N = 3.5;
a = 0.28; b = 9.61; km = 1; kM = 10;
theta_au = asin(H / d_au);
P_los_au = 1 / (1 + a * exp(-b * (theta_au - a)));
alpha_au = (alpha_L - alpha_N) * P_los_au + alpha_N;
L_au     = d_au^(-alpha_au);
K_au     = db2pow(km + (kM - km) * 2 * theta_au / pi);

N = 1e6;
lambda_au = sqrt(2 * K_au);

% Generate and normalize
X_raw = ncx2rnd(2, lambda_au, [1, N]);
X     = (X_raw / mean(X_raw)) * L_au;

fprintf('Target E[X] = L_au = %.6f\n', L_au);
fprintf('Simulated   mean(X) = %.6f\n', mean(X));
fprintf('Relative error = %.4f%%\n\n', abs(mean(X) - L_au) / L_au * 100);

% Rayleigh (G2G)
L_ae = d_au^(-alpha_N);   % G2G path loss
V    = L_ae * exprnd(1, [1, N]);
fprintf('Target E[V] = L_ae = %.6f\n', L_ae);
fprintf('Simulated   mean(V) = %.6f\n', mean(V));
fprintf('Relative error = %.4f%%\n', abs(mean(V) - L_ae) / L_ae * 100);

assert(abs(mean(X) - L_au) / L_au < 0.01, 'Channel mean error > 1%');
assert(abs(mean(V) - L_ae) / L_ae < 0.01, 'Rayleigh channel mean error > 1%');
fprintf('\nAll assertions passed.\n');

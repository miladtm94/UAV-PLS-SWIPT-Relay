%% test_marcumq_approximation.m
%
%  Validates the series expansion approximation of the first-order
%  Marcum-Q function Q_1(a, b(x)) that appears in the analytical
%  connection probability expression (see Sec. III of the paper).
%
%  In the CP derivation, the integral
%
%      CP = integral_0^inf  Q_1(sqrt(2*K_au), sqrt(2*(1+K_au)*(A + B/x)))
%                           * f_Y(x)  dx
%
%  requires evaluating Q_1 at a vector of b-values parameterised by x.
%  This script compares:
%       (a) The R-term series expansion (out1) used in the analytical formula
%       (b) MATLAB's built-in marcumq(a, b) (out2)
%
%  A relative error table is printed; values below 1% confirm the
%  expansion is numerically accurate for the given parameter set.
%
%  Reference: MarqumQfunc_Approximation_Test.m (original debug script)
%
%  Usage:
%       test_marcumq_approximation
%
% -------------------------------------------------------------------------

clc; clear;

fprintf('=== Marcum-Q Function Series Approximation Test ===\n\n');

%% System parameters (matching the main simulation setup)
D      = 5;   H = 1;
a_itu  = 0.28; b_itu = 9.61;  % ITU-R urban parameters
km     = 1;   kM = 8;
alpha_L = 2;  alpha_N = 3.5;
N0     = 1e-2; eff = 0.7; zeta = 2;
beta   = 0.7;
Rt     = 1;   delta_t = 2^(2*Rt) - 1;
Pa     = db2pow(10);          % Alice transmit power (10 dBW)

%% Node positions
t_A = [0,     0, 0];
t_U = [D/5,   0, H];
t_B = [D,     0, 0];
t_E = [4*D/5, 1, 0];

%% Distance and channel statistics
d_au = norm(t_A - t_U);
d_ub = norm(t_U - t_B);

% --- U-B link (Rician) ---
theta_ub    = asin(H / d_ub);
P_los_ub    = 1 / (1 + a_itu * exp(-b_itu*(theta_ub - a_itu)));
alpha_ub    = (alpha_L - alpha_N)*P_los_ub + alpha_N;
L_ub        = d_ub^(-alpha_ub);
K_ub        = db2pow(km + (kM - km)*2*theta_ub/pi);

% --- A-U link (Rician) ---
theta_au    = asin(H / d_au);
P_los_au    = 1 / (1 + a_itu * exp(-b_itu*(theta_au - a_itu)));
alpha_au    = (alpha_L - alpha_N)*P_los_au + alpha_N;
L_au        = d_au^(-alpha_au);
K_au        = db2pow(km + (kM - km)*2*theta_au/pi);

%% CP threshold constants (see Eq. in paper)
A = (1 - beta + zeta) * N0 * delta_t / (1 - beta) / Pa / L_au;
B = N0 * delta_t / eff / beta / Pa / L_au / L_ub;

fprintf('  Threshold constants:  A = %.4e,  B = %.4e\n\n', A, B);

%% Marcum-Q arguments
a1 = sqrt(2 * K_au);                          % fixed argument
x  = (0.1 : 0.1 : 5.0)';                      % relay-channel variable
a2 = sqrt(2*(1 + K_au)*(A + B./x));           % varying argument

%% --- Series Expansion (R-term) ---
R    = 50;
out1 = zeros(size(x));
for r = 0:R
    for u = 0:r
        out1 = out1 + ...
            gamma(R+r) .* R^(1-2*r) .* a1^(2*r) .* a2.^(2*u) ...
            .* exp(-0.5*(a1^2 + a2.^2)) ...
            ./ factorial(u) ./ factorial(r) ./ gamma(R-r+1) ./ 2^(r+u);
    end
end

%% --- MATLAB built-in ---
out2 = marcumq(a1, a2);

%% --- Relative error ---
rel_err = abs(out1 - out2) ./ (out2 + eps);

%% --- Report ---
fprintf('  x       Series(R=50)  MATLAB marcumq  RelErr(%%)\n');
fprintf('  ------  ------------  --------------  ---------\n');
for k = 1:5:length(x)
    fprintf('  %5.2f   %10.6f    %10.6f      %.4f\n', ...
        x(k), out1(k), out2(k), rel_err(k)*100);
end

max_err = max(rel_err) * 100;
fprintf('\n  Max relative error over %d sample points: %.4f %%\n', ...
    length(x), max_err);

if max_err < 1.0
    fprintf('  PASS: Series approximation error < 1%% threshold.\n');
else
    fprintf('  WARN: Max error %.2f%% exceeds 1%% threshold. Increase R.\n', max_err);
end

%% --- Optional: plot ---
figure;
plot(x, out2, '-b', 'LineWidth', 2, 'DisplayName', 'MATLAB marcumq'); hold on;
plot(x, out1, '--r', 'LineWidth', 1.5, 'DisplayName', sprintf('Series (R=%d)', R));
legend('Location','best','FontSize',10,'FontName','Times New Roman');
xlabel('x (relay channel variable)', 'FontSize', 12, 'FontName', 'Times New Roman');
ylabel('Q_1(a, b(x))', 'FontSize', 12, 'FontName', 'Times New Roman');
title('Marcum-Q Approximation Accuracy', 'FontSize', 12, 'FontName', 'Times New Roman');
grid on;

fprintf('\n  Figure saved if output directory is writable.\n');

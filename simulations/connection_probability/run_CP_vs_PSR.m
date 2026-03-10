%% CONNECTION PROBABILITY vs. POWER SPLITTING RATIO (PSR)
%
%  Reproduces the Connection Probability (CP) vs PSR figure showing how
%  the choice of power splitting ratio beta at the UAV affects CP.
%
%  Compares: Monte Carlo simulation, exact integral, and series approximation.
%
%  Reference: Fig. X in the published paper.
%
% -------------------------------------------------------------------------

clc; clear; close all;

addpath(fullfile(fileparts(mfilename('fullpath')), '..', '..', 'src', 'analytics'));

%% ========== Network Topology ============================================
H    = 50;        % UAV altitude (m)
r_au = 40;        % Horizontal distance: Alice to UAV (m)
r_ub = 60;        % Horizontal distance: UAV to Bob (m)
d_au = sqrt(r_au^2 + H^2);   % 3D distance Alice-UAV
d_ub = sqrt(r_ub^2 + H^2);   % 3D distance UAV-Bob

%% ========== Channel Parameters ==========================================
% Non-centrality parameter for A2G links (empirical model: -17.4*log10(H)+29.6)
lambda_au = -17.4 * log10(H) + 29.6;
lambda_ub = -17.4 * log10(H) + 29.6;

% Path-loss models
a_los  = 2.0;   % LoS path-loss exponent
a_nlos = 3.0;   % NLoS path-loss exponent
eta_los_dB  = 1;   eta_los  = db2pow(eta_los_dB);
eta_nlos_dB = 20;  eta_nlos = db2pow(eta_nlos_dB);

% A2G LoS probability model
a_prob = 9.61;  b_prob = 0.28;
P_los_au  = 1 / (1 + a_prob * exp(-b_prob * (asin(H / d_au) - a_prob)));
P_los_ub  = 1 / (1 + a_prob * exp(-b_prob * (asin(H / d_ub) - a_prob)));
L_au = 1 / (eta_los * d_au^a_los);
L_los_ub  = 1 / (eta_los  * d_ub^a_los);
L_nlos_ub = 1 / (eta_nlos * d_ub^a_nlos);
L_ub = P_los_ub * L_los_ub + (1 - P_los_ub) * L_nlos_ub;

%% ========== System Parameters ===========================================
Pa_dBW = 10;                          % Transmit power [dBW]
N0_dBm = -80; N0 = 1e-3 * db2pow(N0_dBm);
Pa  = db2pow(Pa_dBW);
rho = Pa / N0;
e   = 0.5;   % Energy harvesting efficiency
Rt  = 2;
delta_t = 2^(2 * Rt) - 1;

%% ========== Sweep Setup =================================================
n = 30;
N = 1e5;
alphaVec = linspace(1e-8, 1 - 1e-8, n);   % PSR = alpha (beta in other scripts)

Pc_Simulation   = zeros(1, n);
Pc_Exact        = zeros(1, n);
Pc_Approximate  = zeros(1, n);

%% ========== Simulation ==================================================
X = ncx2rnd(2, lambda_au, [1, N]);
Y = ncx2rnd(2, lambda_ub, [1, N]);

for i = 1:n
    alpha = alphaVec(i);
    Gamma_AB = (e * alpha * (1 - alpha) * rho^2 * L_au^2 * X.^2 * L_ub .* Y) ./ ...
               (e * alpha * rho * L_au * L_ub .* X .* Y + (1 - alpha) * rho * L_au .* X + 1);
    Ct = 0.5 * log2(1 + Gamma_AB);
    Pc_Simulation(i) = mean(Ct > Rt);
end

%% ========== Exact Analytical ============================================
xmin = 0; xmax = 20;
for i = 1:n
    alpha = alphaVec(i);
    A = delta_t / (1 - alpha) / rho / L_au;
    B = delta_t / rho / L_au / e / alpha / L_ub;
    Func = @(x) 0.5 .* exp(-(x + lambda_ub) / 2) .* besseli(0, sqrt(lambda_ub .* x)) ...
                    .* (1 - marcumq(sqrt(lambda_au), sqrt(A + B ./ x)));
    Pc_Exact(i) = 1 - integral(Func, xmin, xmax);
end

%% ========== Closed-Form Approximation ===================================
D_ord = 1;
xmax  = 20;
for i = 1:n
    alpha = alphaVec(i);
    A = delta_t / (1 - alpha) / rho / L_au;
    B = delta_t / rho / L_au / e / alpha / L_ub;
    R_ord = ceil(lambda_ub * xmax / 2);
    Func = 0;
    for d = 0:D_ord
        for u = 0:d
            for v = 0:u
                for r = 0:R_ord
                    Func = Func + ...
                        (gamma(d + D_ord) * D_ord^(1 - 2*d) * lambda_au^d * factorial(d) ...
                         * lambda_ub^r * A^(u - v) * R_ord^(1 - 2*r) * nchoosek(u, v) ...
                         * gamma(R_ord + r) * 2^(1 - u - 2*r - d)) ...
                        / ((gamma(d + 1))^2 * gamma(D_ord - d + 1) * factorial(u) ...
                           * (gamma(r + 1))^2 * gamma(R_ord - r + 1)) ...
                        * B^((r + v + 1) / 2) * besselk(r - v + 1, sqrt(B));
                end
            end
        end
    end
    Pc_Approximate(i) = 0.5 * exp(-(lambda_ub + lambda_au + A) / 2) * Func;
end

%% ========== Plot ========================================================
figure('Name', 'Connection Probability vs PSR');
plot(alphaVec, Pc_Exact,      '-or', 'LineWidth', 1.5, 'DisplayName', 'Exact (Integral)');
hold on;
plot(alphaVec, Pc_Approximate,'--sb','LineWidth', 1.5, 'DisplayName', 'Approximate (Series)');
plot(alphaVec, Pc_Simulation, '*k',  'MarkerSize', 6,  'DisplayName', 'Monte Carlo');

xlabel('Power Splitting Ratio (\alpha)',  'FontSize', 14, 'FontName', 'Times New Roman');
ylabel('Connection Probability',           'FontSize', 14, 'FontName', 'Times New Roman');
legend('show', 'Location', 'best', 'FontSize', 11);
grid on; box on;

%% ========== Save ========================================================
results_dir = fullfile(fileparts(mfilename('fullpath')), '..', '..', 'results', 'figures');
if ~exist(results_dir, 'dir'); mkdir(results_dir); end
savefig(fullfile(results_dir, 'CP_vs_PSR.fig'));
print( fullfile(results_dir, 'CP_vs_PSR.png'), '-dpng', '-r300');

MSE = mean(abs(Pc_Simulation - Pc_Approximate) ./ Pc_Simulation);
fprintf('Mean Absolute Relative Error (Approximate vs Simulation): %.3f%%\n', MSE * 100);

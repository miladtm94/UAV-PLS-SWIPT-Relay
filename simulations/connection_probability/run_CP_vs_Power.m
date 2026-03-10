%% CONNECTION PROBABILITY vs. NETWORK TRANSMIT POWER
%
%  Reproduces the Connection Probability (CP) figure comparing:
%    - Monte Carlo simulation (exact)
%    - Closed-form integral expression (exact analytical)
%    - Closed-form series approximation (analytical approximation)
%
%  System: Alice (A) --> UAV relay (U) --> Bob (B)
%  Channel model: Rician fading for A2G links (elevation-dependent K-factor),
%                 Rayleigh fading for G2G links.
%  SWIPT: Power Splitting Receiver (PSR) at UAV with ratio beta.
%  Power Allocation: Equal Power Allocation (lambda = 0.5).
%
%  Reference: Fig. X in the published paper (IEEE Access, 2024).
%
%  Run time: ~2-5 minutes (N = 1e5 Monte Carlo samples).
%
%  Dependencies: src/analytics/coeff.m, MATLAB Statistics Toolbox.
%
%  Output: Figure saved to results/figures/CP_vs_Power.fig and .png
%
% -------------------------------------------------------------------------
% Author: [Author Name(s)]
% Institution: [Institution]
% Paper DOI: 10.1109/ACCESS.XXXXXXX
% -------------------------------------------------------------------------

clc; clear; close all;

%% ========== Add source paths ===========================================
addpath(fullfile(fileparts(mfilename('fullpath')), '..', '..', 'src', 'analytics'));
addpath(fullfile(fileparts(mfilename('fullpath')), '..', '..', 'src', 'channel'));

%% ========== Network Topology ============================================
D  = 5;     % Alice-to-Bob horizontal distance (normalized units)
H  = 1;     % UAV altitude (normalized units)

% Node positions [x, y, z]
t_A = [0,   0, 0];      % Alice (source)
t_U = [D/4, 0, H];      % UAV relay
t_B = [D,   0, 0];      % Bob (destination)

% Inter-node distances
d_au = norm(t_A - t_U); % Alice  -> UAV
d_ub = norm(t_U - t_B); % UAV   -> Bob

%% ========== Channel Parameters ==========================================
alpha_L = 2;    % LoS path-loss exponent
alpha_N = 3.5;  % NLoS path-loss exponent

% Urban LoS probability model parameters
a = 0.28;   % [ITU-R / empirical]
b = 9.61;

% Rician K-factor range [dB]
km = 1;    % Minimum K (near-horizontal links)
kM = 10;   % Maximum K (near-vertical links)

% Compute A2G channel statistics (elevation-dependent)
theta_au    = asin(H / d_au);
P_los_au    = 1 / (1 + a * exp(-b * (theta_au - a)));
alpha_au    = (alpha_L - alpha_N) * P_los_au + alpha_N;
L_au        = d_au^(-alpha_au);              % Large-scale path loss
K_au        = db2pow(km + (kM - km) * 2 * theta_au / pi);  % K-factor

theta_ub    = asin(H / d_ub);
P_los_ub    = 1 / (1 + a * exp(-b * (theta_ub - a)));
alpha_ub    = (alpha_L - alpha_N) * P_los_ub + alpha_N;
L_ub        = d_ub^(-alpha_ub);
K_ub        = db2pow(km + (kM - km) * 2 * theta_ub / pi);

lambda_au   = sqrt(2 * K_au);  % ncx2 non-centrality parameter
lambda_ub   = sqrt(2 * K_ub);

%% ========== SWIPT & System Parameters ==================================
N0   = 1e-2;   % Noise power (linear)
eta  = 0.5;    % Energy harvesting efficiency
zeta = 2;      % Pilot power to noise ratio (Np/N0)
beta = 0.5;    % Power splitting ratio (0 < beta < 1)
Rt   = 2;      % Target transmission rate [bits/s/Hz]
delta_t = 2^(2 * Rt) - 1;  % SNR threshold

%% ========== Simulation Setup ============================================
n  = 30;    % Number of power levels
N  = 1e5;   % Monte Carlo samples per point

P_dBW  = linspace(-5, 30, n);
Pa_Vec = db2pow(P_dBW) / 2;   % Alice's power (EPA: lambda = 0.5)
Pb_Vec = db2pow(P_dBW) / 2;   % Bob's jamming power

Pc_Simulation   = zeros(1, n);
Pc_Exact        = zeros(1, n);
Pc_Approximate  = zeros(1, n);

%% ========== Monte Carlo Simulation =====================================
fprintf('Running Monte Carlo simulation...\n');
for i = 1:n
    X = ncx2rnd(2, lambda_au, [1, N]);
    Y = ncx2rnd(2, lambda_ub, [1, N]);

    Pa = Pa_Vec(i);
    Pb = Pb_Vec(i);

    % Effective noise due to pilot contamination
    epsilon = zeta * N0^2 ./ (Pa * X + Pb * Y);

    % End-to-end SNR at Bob (SWIPT-based DF relay)
    Gamma_AB = (eta * beta * (1 - beta) * Pa .* L_au .* X .* L_ub .* Y) ./ ...
               (eta * beta * (1 - beta + zeta) * L_ub .* Y * N0 + (1 - beta) * N0 + epsilon);

    Ct = 0.5 * log2(1 + Gamma_AB);
    Pc_Simulation(i) = mean(Ct > Rt);
end

%% ========== Exact Analytical (Numerical Integration) ===================
fprintf('Computing exact analytical (integral form)...\n');
xmax = 20;   % Integration upper limit (sufficient for ncx2 CDF)
for i = 1:n
    Pa = Pa_Vec(i);
    Pb = Pb_Vec(i);
    A  = (1 - beta + zeta) * delta_t / (1 - beta) / (Pa / N0) / L_au;
    B  = delta_t / (Pa / N0) / L_au / eta / beta / L_ub;

    Func = @(x) 0.5 .* exp(-(x + lambda_ub) / 2) .* besseli(0, sqrt(lambda_ub .* x)) ...
                    .* (1 - marcumq(sqrt(lambda_au), sqrt(A + B ./ x)));
    Pc_Exact(i) = 1 - integral(Func, 0, xmax);
end

%% ========== Closed-Form Series Approximation ============================
fprintf('Computing closed-form series approximation...\n');
D_ord = 15;  % Series truncation order for ncx2 expansion
R_ord = 15;  % Series truncation order for Rician expansion

for i = 1:n
    Pa = Pa_Vec(i);
    Pb = Pb_Vec(i);
    A  = (1 - beta + zeta) * delta_t / (1 - beta) / (Pa / N0) / L_au;
    B  = delta_t / (Pa / N0) / L_au / eta / beta / L_ub;

    Func = 0;
    for d = 0:D_ord
        for u = 0:d
            for s = 0:u
                for r = 0:R_ord
                    Func = Func + 2 * ( ...
                        (gamma(d + D_ord) * D_ord^(1 - 2*d) * K_au^d * (1 + K_au)^u * ...
                         B^(u - s) * A^s * R_ord^(1 - 2*r) * nchoosek(u, s) * gamma(R_ord + r)) ...
                        / (gamma(d + 1) * gamma(D_ord - d + 1) * gamma(u - s + 1) * ...
                           gamma(s + 1) * (gamma(r + 1))^2 * gamma(R_ord - r + 1)) ...
                        * (K_ub * (1 + K_ub))^r ...
                        * ((1 + K_au) * B / (1 + K_ub))^((r - (u - s) + 1) / 2) ...
                        * besselk(r - (u - s) + 1, 2 * sqrt((1 + K_au) * B * (1 + K_ub))) );
                end
            end
        end
    end
    Pc_Approximate(i) = (1 + K_ub) * exp(-K_au - K_ub - (1 + K_au) * A) * Func;
end

%% ========== Plot Results ================================================
figure('Name', 'Connection Probability vs Transmit Power', 'NumberTitle', 'off');
plot(P_dBW, Pc_Exact,       '-or',  'LineWidth', 1.5, 'DisplayName', 'Exact (Integral)');
hold on;
plot(P_dBW, Pc_Approximate, '--sb', 'LineWidth', 1.5, 'DisplayName', 'Approximate (Series)');
plot(P_dBW, Pc_Simulation,  '*k',   'MarkerSize', 6,  'DisplayName', 'Monte Carlo Simulation');

xlabel('Network Transmit Power P [dBW]', 'FontSize', 14, 'FontName', 'Times New Roman');
ylabel('Connection Probability',         'FontSize', 14, 'FontName', 'Times New Roman');
legend('show', 'Location', 'best', 'FontSize', 11);
grid on; box on; axis tight;
set(gca, 'FontSize', 12);

%% ========== Save Results ================================================
results_dir = fullfile(fileparts(mfilename('fullpath')), '..', '..', 'results', 'figures');
if ~exist(results_dir, 'dir'); mkdir(results_dir); end

savefig(fullfile(results_dir, 'CP_vs_Power.fig'));
print( fullfile(results_dir, 'CP_vs_Power.png'), '-dpng', '-r300');

% Save numerical data
data_dir = fullfile(fileparts(mfilename('fullpath')), '..', '..', 'data', 'processed');
if ~exist(data_dir, 'dir'); mkdir(data_dir); end

save(fullfile(data_dir, 'CP_vs_Power.mat'), ...
     'P_dBW', 'Pc_Simulation', 'Pc_Exact', 'Pc_Approximate', ...
     'D', 'H', 'alpha_L', 'alpha_N', 'beta', 'eta', 'zeta', 'Rt', 'N');

fprintf('\nResults saved to results/figures/ and data/processed/\n');

%% ========== Accuracy Report =============================================
MSE = mean(Pc_Exact .* (Pc_Exact - Pc_Approximate).^2) / mean(Pc_Exact);
fprintf('Weighted MSE (Exact vs Approximate): %.5f%%\n', MSE * 100);
fprintf('Table: [P_dBW, Pc_Exact, Pc_Approximate, Pc_Simulation]\n');
disp([P_dBW', Pc_Exact', Pc_Approximate', Pc_Simulation']);

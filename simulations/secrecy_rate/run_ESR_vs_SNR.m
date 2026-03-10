%% ERGODIC SECRECY RATE vs. NETWORK TRANSMIT SNR
%
%  Reproduces the Average Secrecy Rate (ASR) vs transmit SNR figure,
%  comparing the analytical lower bound (closed-form) against Monte Carlo
%  simulation.
%
%  System:  Alice (A) --[Rician]--> UAV relay (U) --[Rician]--> Bob (B)
%           Eavesdropper (E) receives from both A (Rayleigh) and U (Rician).
%  Security mechanism: Bob transmits cooperative jamming signal.
%  SWIPT at UAV: Power Splitting Receiver (PSR) with ratio beta.
%
%  The ESR lower bound uses the Jensen-type approximation:
%       Rs_LB = 0.5 * [log(1 + exp(T1)) - log(1 + T2)]
%  where T1 approximates E[log Gamma_AB] and T2 approximates E[Gamma_E].
%
%  Reference: Fig. X in the published paper (IEEE Access, 2024).
%
% -------------------------------------------------------------------------

clc; clear; close all;

addpath(fullfile(fileparts(mfilename('fullpath')), '..', '..', 'src', 'analytics'));

%% ========== Network Topology ============================================
d  = 2;    H = 1;
Phi_E = pi/2;   r_E = d/2;

d_au = sqrt((d/2)^2 + H^2);    % Alice  -> UAV (symmetric placement)
d_bu = sqrt((d/2)^2 + H^2);    % Bob   -> UAV
d_ue = sqrt(r_E^2  + H^2);     % UAV   -> Eve
d_ae = sqrt((d/2)^2 + r_E^2 - r_E * d * cos(Phi_E));       % Alice -> Eve
d_be = sqrt((d/2)^2 + r_E^2 - r_E * d * cos(pi - Phi_E));  % Bob   -> Eve

%% ========== Channel Parameters ==========================================
alpha_G = 3.5;  % G2G path-loss exponent
alpha_L = 2.0;  % LoS path-loss exponent
alpha_N = 3.0;  % NLoS path-loss exponent
eta = 0.7;      % Energy harvesting efficiency
beta = 0.5;     % Power splitting ratio

% Rician channel non-centrality parameters
mu1_x = d_au^(-alpha_L);  mu2_x = d_au^(-alpha_N);
mu1_y = d_bu^(-alpha_L);  mu2_y = d_bu^(-alpha_N);
mu1_z = d_ue^(-alpha_L);  mu2_z = d_ue^(-alpha_N);
lambda_x = mu1_x^2 + mu2_x^2;
lambda_y = mu1_y^2 + mu2_y^2;
lambda_z = mu1_z^2 + mu2_z^2;

% Rayleigh channel means (G2G: Eve links)
mu2_v = d_ae^(-alpha_G);   % Alice->Eve mean power
mu2_w = d_be^(-alpha_G);   % Bob->Eve mean power

%% ========== System Parameters ===========================================
N0_dBm = -20;  N0 = 1e-3 * db2pow(N0_dBm);
n  = 30;   N  = 1e4;
P_dBW = linspace(0, 20, n);
rhoVec = db2pow(P_dBW) / N0;

%% ========== Channel Simulation ==========================================
X = ncx2rnd(2, lambda_x, [n, N]);
Y = ncx2rnd(2, lambda_y, [n, N]);
Z = ncx2rnd(2, lambda_z, [n, N]);
V = chi2rnd(mu2_v, [n, N]);   % Rayleigh power ~ scaled chi2(1)
W = chi2rnd(mu2_w, [n, N]);

%% ========== Monte Carlo Simulation =====================================
fprintf('Running Monte Carlo simulation...\n');
RhoVec  = repmat(rhoVec', 1, N);

Gamma_AB = ((1 - beta) .* RhoVec .* X .* Y) ./ (Y + (1 - beta) / eta / beta);
Gamma_e1 = RhoVec .* V ./ (RhoVec .* W + 1);
Gamma_e2 = (1 - beta) .* RhoVec .* X .* Z ./ ...
           (((1 - beta) .* RhoVec .* Y + 1) .* Z + (1 - beta) / eta / beta);
Gamma_E  = Gamma_e1 + Gamma_e2;

Rsec    = 0.5 * log2((1 + Gamma_AB) ./ (1 + Gamma_E));
ESR_sim = mean(Rsec, 2);

%% ========== Analytical Lower Bound =====================================
fprintf('Computing analytical lower bound...\n');
mx = lambda_x + 2;   % E[X] for normalized ncx2(2,lambda_x)
my = lambda_y + 2;
mz = lambda_z + 2;
mw = mu2_w;   mv = mu2_v;

R_terms = 10;   % Series truncation for g1, g2 approximations

ESR_analytical = zeros(1, n);
for i = 1:n
    rho = rhoVec(i);

    % T1: lower bound on E[log Gamma_AB]
    T1 = log((1 - beta) * rho) + g1(lambda_x, R_terms) + g1(lambda_y, R_terms) ...
         - g2(lambda_y, (1 - beta) / eta / beta, R_terms);

    % T2: upper bound on E[Gamma_E]
    L2 = rho * (1 - beta) * mx / ...
         ((1 - beta) * rho * my + ...
          2 * (1 - beta) * exp(-lambda_z/4) * whittakerM(-0.5, 0, lambda_z/2) ...
          / eta / beta / sqrt(2 * lambda_z) + 1);
    mww = mw * rho;
    L1  = rho * mv * (1 / mww) * exp(1 / mww) * expint(1 / mww);
    T2  = L1 + L2;

    ESR_analytical(i) = max((1/2/log(2)) * (log(1 + exp(T1)) - log(1 + T2)), 0);
end

%% ========== Plot ========================================================
figure('Name', 'Ergodic Secrecy Rate vs Transmit SNR');
plot(P_dBW, ESR_analytical, '-or', 'LineWidth', 1.5, 'DisplayName', 'Analytical LB');
hold on;
plot(P_dBW, ESR_sim,        '--b', 'LineWidth', 1.5, 'DisplayName', 'Simulation (Exact)');

xlabel('Network Transmit Power P [dBW]', 'FontSize', 14, 'FontName', 'Times New Roman');
ylabel('Ergodic Secrecy Rate [bits/s/Hz]', 'FontSize', 14, 'FontName', 'Times New Roman');
legend('show', 'Location', 'best', 'FontSize', 11);
grid on; box on;

%% ========== Save ========================================================
results_dir = fullfile(fileparts(mfilename('fullpath')), '..', '..', 'results', 'figures');
if ~exist(results_dir, 'dir'); mkdir(results_dir); end
savefig(fullfile(results_dir, 'ESR_vs_SNR.fig'));
print( fullfile(results_dir, 'ESR_vs_SNR.png'), '-dpng', '-r300');

data_dir = fullfile(fileparts(mfilename('fullpath')), '..', '..', 'data', 'processed');
if ~exist(data_dir, 'dir'); mkdir(data_dir); end
save(fullfile(data_dir, 'ESR_vs_SNR.mat'), ...
     'P_dBW', 'ESR_sim', 'ESR_analytical', 'd', 'H', 'beta', 'eta', 'N');

MSE = 100 * mean(abs(ESR_analytical - ESR_sim').^2);
fprintf('MSE (Analytical vs Simulation): %.5f%%\n', MSE);

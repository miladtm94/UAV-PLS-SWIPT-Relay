%% ERGODIC SECRECY RATE vs. POWER SPLITTING RATIO (PSR) — Main Result
%
%  This is the principal simulation script of the paper. It sweeps the
%  Power Splitting Ratio (beta) and computes all four key performance
%  metrics:
%    1. Connection Probability (CP)
%    2. Intercept Probability (IP)
%    3. Ergodic Secrecy Rate (ESR)
%    4. Secrecy Outage Probability (SOP)
%
%  Additionally, the analytical lower bound on the ESR is compared against
%  Monte Carlo simulation, validating the closed-form expression derived
%  in the paper.
%
%  Reference: Fig. X (subplot figure) and Fig. Y (ESR vs beta) in the paper.
%
% -------------------------------------------------------------------------

clc; clear; close all;

addpath(fullfile(fileparts(mfilename('fullpath')), '..', '..', 'src', 'analytics'));

%% ========== Network Topology ============================================
D = 10; H = 1.5;
t_A = [0,     0, 0];
t_U = [D/5,   0, H];
t_R = [D/5,   0, 0];   % Ground relay reference point
t_B = [D,     0, 0];
t_E = [4*D/5, 1, 0];

d_au = norm(t_A - t_U);  d_bu = norm(t_U - t_B);
d_ar = norm(t_A - t_R);  d_rb = norm(t_R - t_B);
d_ae = norm(t_A - t_E);  d_be = norm(t_B - t_E);
d_ue = norm(t_U - t_E);  d_re = norm(t_R - t_E);

%% ========== Channel Parameters ==========================================
alpha_G = 3.5;  alpha_L = 2; alpha_N = 3.5;
eta = 0.5;      Pa = 10;     Pb = 10;
N0 = 1e-2;      beta0 = 1;

% Rician non-centrality parameters
mu1_x = beta0 * d_au^(-alpha_L);  mu2_x = beta0 * d_au^(-alpha_N);
mu1_y = beta0 * d_bu^(-alpha_L);  mu2_y = beta0 * d_bu^(-alpha_N);
mu1_z = beta0 * d_ue^(-alpha_L);  mu2_z = beta0 * d_ue^(-alpha_N);
lambda_x = mu1_x^2 + mu2_x^2;
lambda_y = mu1_y^2 + mu2_y^2;
lambda_z = mu1_z^2 + mu2_z^2;

% Rayleigh (G2G) mean powers
mu2_v = beta0 * d_ae^(-alpha_G);
mu2_w = beta0 * d_be^(-alpha_G);

%% ========== Target Rates ================================================
R_t = 5;    % Target transmission rate [bits/s/Hz]
R_s = 3;    % Target secrecy rate
delta_t = 2^(2*R_t) - 1;
delta_e = 2^(2*(R_t - R_s)) - 1;

%% ========== Simulation Setup ============================================
n  = 20;
N  = 1e5;
betaVec  = linspace(0.01, 0.99, n);
BetaVec  = repmat(betaVec', 1, N);

%% ========== Channel Generation ==========================================
X = ncx2rnd(2, lambda_x, [n, N]);
Y = ncx2rnd(2, lambda_y, [n, N]);
Z = ncx2rnd(2, lambda_z, [n, N]);
V = chi2rnd(mu2_v, [n, N]);
W = chi2rnd(mu2_w, [n, N]);

%% ========== Monte Carlo: All Metrics ====================================
fprintf('Computing all metrics via Monte Carlo...\n');
Gamma_AB = ((1 - BetaVec) .* (Pa/N0) .* X .* Y) ./ ...
           (Y + (1 - BetaVec) ./ eta ./ BetaVec);
Gamma_e1 = (Pa/N0) .* V ./ ((Pb/N0) .* W + 1);
Gamma_e2 = (1 - BetaVec) .* (Pa/N0) .* X .* Z ./ ...
           (((1 - BetaVec) .* (Pb/N0) .* Y + 1) .* Z + (1 - BetaVec) ./ eta ./ BetaVec);
Gamma_E  = Gamma_e1 + Gamma_e2;

Rsec    = max(0.5 * log2((1 + Gamma_AB) ./ (1 + Gamma_E)), 0);
CP_sim  = mean(Gamma_AB > delta_t, 2);
IP_sim  = mean(Gamma_E  > delta_e, 2);
ESR_sim = mean(Rsec, 2);
SOP_sim = mean(Rsec < R_s, 2);

%% ========== Analytical ESR Lower Bound ==================================
fprintf('Computing analytical ESR lower bound...\n');
mx = lambda_x + 2;  my = lambda_y + 2;  mz = lambda_z + 2;
mw = mu2_w;         mv = mu2_v;

R_terms = 10;
ESR_analytical = zeros(n, 1);
for i = 1:n
    beta = betaVec(i);
    T1 = log((1 - beta) * Pa / N0) ...
         + g1(lambda_x, R_terms) + g1(lambda_y, R_terms) ...
         - g2(lambda_y, (1 - beta) / eta / beta, R_terms);
    L2 = (Pa/N0) * (1 - beta) * mx / ...
         ((1 - beta) * (Pb/N0) * my + ...
          2*(1-beta)*exp(-lambda_z/4)*whittakerM(-0.5,0,lambda_z/2) ...
          / eta/beta/sqrt(2*lambda_z) + 1);
    mww = mw * (Pb/N0);
    L1  = (Pa/N0) * mv * (1/mww) * exp(1/mww) * expint(1/mww);
    T2  = L1 + L2;
    ESR_analytical(i) = max((1/2/log(2)) * (log(1 + exp(T1)) - log(1 + T2)), 0);
end

%% ========== Plot 1: All Metrics vs PSR ==================================
figure('Name', 'All Performance Metrics vs PSR', 'Position', [100 100 900 700]);

subplot(2, 2, 1);
plot(betaVec, CP_sim, '-r', 'LineWidth', 1.5);
xlabel('PSR (\beta)', 'FontName', 'Times New Roman', 'FontSize', 12);
ylabel('Connection Probability', 'FontName', 'Times New Roman', 'FontSize', 12);
grid on;

subplot(2, 2, 2);
plot(betaVec, IP_sim, '-b', 'LineWidth', 1.5);
xlabel('PSR (\beta)', 'FontName', 'Times New Roman', 'FontSize', 12);
ylabel('Intercept Probability', 'FontName', 'Times New Roman', 'FontSize', 12);
grid on;

subplot(2, 2, 3);
plot(betaVec, ESR_sim, '-g', 'LineWidth', 1.5);
xlabel('PSR (\beta)', 'FontName', 'Times New Roman', 'FontSize', 12);
ylabel('Ergodic Secrecy Rate [bits/s/Hz]', 'FontName', 'Times New Roman', 'FontSize', 12);
grid on;

subplot(2, 2, 4);
plot(betaVec, SOP_sim, '-k', 'LineWidth', 1.5);
xlabel('PSR (\beta)', 'FontName', 'Times New Roman', 'FontSize', 12);
ylabel('Secrecy Outage Probability', 'FontName', 'Times New Roman', 'FontSize', 12);
grid on;

sgtitle('Performance Metrics vs Power Splitting Ratio', ...
        'FontName', 'Times New Roman', 'FontSize', 14);

%% ========== Plot 2: ESR Analytical vs Simulation ========================
figure('Name', 'ESR: Analytical Lower Bound vs Simulation');
plot(betaVec, ESR_analytical, '-or', 'LineWidth', 1.5, 'DisplayName', 'Analytical LB');
hold on;
plot(betaVec, ESR_sim,        '--b', 'LineWidth', 1.5, 'DisplayName', 'Simulation');

xlabel('Power Splitting Ratio (\beta)',    'FontSize', 14, 'FontName', 'Times New Roman');
ylabel('Ergodic Secrecy Rate [bits/s/Hz]', 'FontSize', 14, 'FontName', 'Times New Roman');
legend('show', 'Location', 'best', 'FontSize', 11);
grid on; box on;

MSE = 100 * mean(abs(ESR_analytical - ESR_sim).^2);
fprintf('MSE (Analytical vs Simulation): %.5f%%\n', MSE);

%% ========== Save ========================================================
results_dir = fullfile(fileparts(mfilename('fullpath')), '..', '..', 'results', 'figures');
if ~exist(results_dir, 'dir'); mkdir(results_dir); end
savefig(1, fullfile(results_dir, 'All_Metrics_vs_PSR.fig'));
print(1, fullfile(results_dir, 'All_Metrics_vs_PSR.png'), '-dpng', '-r300');
savefig(2, fullfile(results_dir, 'ESR_vs_PSR_Analytical.fig'));
print(2, fullfile(results_dir, 'ESR_vs_PSR_Analytical.png'), '-dpng', '-r300');

data_dir = fullfile(fileparts(mfilename('fullpath')), '..', '..', 'data', 'processed');
if ~exist(data_dir, 'dir'); mkdir(data_dir); end
save(fullfile(data_dir, 'ESR_vs_PSR.mat'), ...
     'betaVec', 'CP_sim', 'IP_sim', 'ESR_sim', 'SOP_sim', 'ESR_analytical', ...
     'D', 'H', 'Pa', 'Pb', 'N0', 'eta', 'R_t', 'R_s', 'N');

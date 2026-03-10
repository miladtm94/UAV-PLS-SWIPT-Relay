%% ERGODIC SECRECY RATE vs. POWER ALLOCATION FACTOR (PAF)
%
%  Demonstrates the effect of the Power Allocation Factor (lambda) on the
%  Average Secrecy Rate. Alice allocates fraction lambda to the information
%  relay link and (1 - lambda) to cooperative jamming via Bob.
%
%  This figure helps identify the optimal lambda* that maximizes the ESR.
%
%  Reference: Fig. X in the published paper.
%
% -------------------------------------------------------------------------

clc; clear; close all;

addpath(fullfile(fileparts(mfilename('fullpath')), '..', '..', 'src', 'analytics'));
addpath(fullfile(fileparts(mfilename('fullpath')), '..', '..', 'src', 'optimization'));

%% ========== Network Topology ============================================
N_mc  = 1e6;   % Large sample for smooth PAF curve
D  = 2;
lambda_geo = 3/10;   % UAV horizontal position ratio
d1 = lambda_geo * D; d2 = (1 - lambda_geo) * D;
H  = 1;
Phi_E = pi/2;  r_E = D/2;

d_au = sqrt(d1^2 + H^2);   d_bu = sqrt(d2^2 + H^2);
d_ue = sqrt(r_E^2 + H^2);
d_ae = sqrt(d1^2 + r_E^2 - 2*r_E*d1*cos(Phi_E));
d_be = sqrt(d2^2 + r_E^2 - 2*r_E*d2*cos(pi - Phi_E));

%% ========== Channel Parameters ==========================================
alpha_G = 4; alpha_L = 2; alpha_N = 3;
eta = 0.5; beta = 0.5;
N0 = 1e-4; beta0 = 1;   % beta0: path-loss reference

% Rician parameters
mu1_x = beta0 * d_au^(-alpha_L);  mu2_x = beta0 * d_au^(-alpha_N);
mu1_y = beta0 * d_bu^(-alpha_L);  mu2_y = beta0 * d_bu^(-alpha_N);
mu1_z = beta0 * d_ue^(-alpha_L);  mu2_z = beta0 * d_ue^(-alpha_N);
lambda_x = mu1_x^2 + mu2_x^2;
lambda_y = mu1_y^2 + mu2_y^2;
lambda_z = mu1_z^2 + mu2_z^2;
mu2_v = beta0 * d_ae^(-alpha_G);
mu2_w = beta0 * d_be^(-alpha_G);

P = 10 * log(20);   % Network transmit power (linear)
rho = P / N0;

%% ========== PAF Sweep ===================================================
n = 30;
lambdaVec = linspace(0, 1, n);
LambdaVec = repmat(lambdaVec', 1, N_mc);
Pa = LambdaVec * P;
Pb = (1 - LambdaVec) * P;

%% ========== Channel Generation ==========================================
X = ncx2rnd(2, lambda_x, [n, N_mc]);
Y = ncx2rnd(2, lambda_y, [n, N_mc]);
Z = ncx2rnd(2, lambda_z, [n, N_mc]);
V = chi2rnd(mu2_v, [n, N_mc]);
W = chi2rnd(mu2_w, [n, N_mc]);

%% ========== SNR Computation =============================================
Gamma_AB = ((1 - beta) .* (Pa/N0) .* X .* Y) ./ (Y + (1 - beta) / eta / beta);
Gamma_e1 = (Pa/N0) .* V ./ ((Pb/N0) .* W + 1);
Gamma_e2 = (1 - beta) .* (Pa/N0) .* X .* Z ./ ...
           (((1 - beta) .* (Pb/N0) .* Y + 1) .* Z + (1 - beta) / eta / beta);
Gamma_E  = Gamma_e1 + Gamma_e2;

Rsec    = max(0.5 * log2((1 + Gamma_AB) ./ (1 + Gamma_E)), 0);
ESR_sim = mean(Rsec, 2);

%% ========== Identify Optimal Lambda =====================================
[ESR_max, idx_opt] = max(ESR_sim);
lambda_opt = lambdaVec(idx_opt);
fprintf('Optimal lambda*  = %.4f\n', lambda_opt);
fprintf('Maximum ESR      = %.4f bits/s/Hz\n', ESR_max);

%% ========== Plot ========================================================
figure('Name', 'Ergodic Secrecy Rate vs Power Allocation Factor');
plot(lambdaVec, ESR_sim, '-or', 'LineWidth', 1.5);
hold on;
xline(lambda_opt, '--k', sprintf('\\lambda^* = %.2f', lambda_opt), ...
      'LabelVerticalAlignment', 'bottom', 'LineWidth', 1);

xlabel('Power Allocation Factor (\lambda)',  'FontSize', 14, 'FontName', 'Times New Roman');
ylabel('Ergodic Secrecy Rate [bits/s/Hz]',   'FontSize', 14, 'FontName', 'Times New Roman');
legend('Simulation', 'Location', 'best', 'FontSize', 11);
grid on; box on;

%% ========== Save ========================================================
results_dir = fullfile(fileparts(mfilename('fullpath')), '..', '..', 'results', 'figures');
if ~exist(results_dir, 'dir'); mkdir(results_dir); end
savefig(fullfile(results_dir, 'ESR_vs_PAF.fig'));
print( fullfile(results_dir, 'ESR_vs_PAF.png'), '-dpng', '-r300');

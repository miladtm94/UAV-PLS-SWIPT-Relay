%% COMPREHENSIVE SIMULATION — All Main Figures
%
%  This master simulation script reproduces the three principal figures of
%  the paper in sequence:
%
%    FIGURE 1: Connection Probability vs Network Transmit Power
%              - UAV Relay (EPA, theory & simulation)
%              - UAV Relay (OPSA — optimal PSR)
%              - Ground Relay (EPA, benchmark)
%
%    FIGURE 2: Intercept Probability vs Network Transmit Power
%              - UAV Relay with Cooperative Jamming (CJ)
%              - UAV Relay without CJ (conventional)
%              - Ground Relay without CJ (benchmark)
%              - Analytical approximation
%
%    FIGURE 3: Average Secrecy Rate vs lambda and beta (3D surface)
%              Joint effect of Power Allocation Factor (lambda) and
%              Power Splitting Ratio (beta) on the ESR.
%
%  Reference: Figs. X, Y, Z in the published paper.
%  Run time: ~20-40 minutes (N = 1e5).
%
% -------------------------------------------------------------------------

clc; clear; close all;

addpath(fullfile(fileparts(mfilename('fullpath')), '..', '..', 'src', 'analytics'));

%% ========== Common System Parameters ====================================
Rt = 0.5; delta_t = 2^(2*Rt) - 1;
Rs = 0.2; Re = Rt - Rs; delta_e = 2^(2*Re) - 1;

D = 10; H = 1.5;
t_A = [0,     0, 0]; t_U = [D/5, 0, H];
t_R = [D/5,   0, 0]; t_B = [D,   0, 0]; t_E = [4*D/5, 1, 0];

d_au = norm(t_A - t_U); d_ub = norm(t_U - t_B);
d_ar = norm(t_A - t_R); d_rb = norm(t_R - t_B);
d_ae = norm(t_A - t_E); d_be = norm(t_B - t_E);
d_ue = norm(t_U - t_E); d_re = norm(t_R - t_E);

alpha_L = 2; alpha_N = 3.5; N0 = 1e-2; eff = 0.7; zeta = 2;
a = 0.28; b = 9.61; km = 1; kM = 10;

theta_ub = asin(H/d_ub); P_los_ub = 1/(1+a*exp(-b*(theta_ub-a)));
L_ub = d_ub^(-(( alpha_L-alpha_N)*P_los_ub+alpha_N));
K_ub = db2pow(km+(kM-km)*2*theta_ub/pi);

theta_au = asin(H/d_au); P_los_au = 1/(1+a*exp(-b*(theta_au-a)));
L_au = d_au^(-((alpha_L-alpha_N)*P_los_au+alpha_N));
K_au = db2pow(km+(kM-km)*2*theta_au/pi);

theta_ue = asin(H/d_ue); P_los_ue = 1/(1+a*exp(-b*(theta_ue-a)));
L_ue = d_ue^(-((alpha_L-alpha_N)*P_los_ue+alpha_N));
K_ue = db2pow(km+(kM-km)*2*theta_ue/pi);

L_ae = d_ae^(-alpha_N); L_be = d_be^(-alpha_N);

n = 20; N = 1e5;
lambda_au = sqrt(2*K_au); lambda_ub = sqrt(2*K_ub); lambda_ue = sqrt(2*K_ue);

X = ncx2rnd(2,lambda_au,[n,N]); X = X./mean(X); X = X*L_au;
Y = ncx2rnd(2,lambda_ub,[n,N]); Y = Y./mean(Y); Y = Y*L_ub;
Z = ncx2rnd(2,lambda_ue,[n,N]); Z = Z./mean(Z); Z = Z*L_ue;
V = L_ae*exprnd(1,[n,N]); W = L_be*exprnd(1,[n,N]);

fprintf('Channels generated.\n');

%% ======================================================================
%% FIGURE 1: Connection Probability vs Transmit Power
%% ======================================================================
fprintf('=== Figure 1: Connection Probability ===\n');
P_dBW  = linspace(5, 30, n);
beta   = 0.5; lambda = 0.5;
PVec   = db2pow(P_dBW);
rhoVec = repmat(PVec', 1, N);
Pa = rhoVec * lambda; Pb = rhoVec * (1 - lambda);

epsilon  = zeta*N0^2 ./ (Pa.*X + Pb.*Y);
Gamma_AB = (eff*beta*(1-beta)*Pa.*X.*Y) ./ ...
           (eff*beta*(1-beta+zeta)*Y*N0 + (1-beta)*N0 + epsilon);
CP_Sim_EPA = mean(0.5*log2(1+Gamma_AB) > Rt, 2);

% Optimal PSR
beta_opt_vec = 1 ./ (1 + sqrt(eff * Y * zeta));
Gamma_AB_opt = (eff*beta_opt_vec.*(1-beta_opt_vec).*Pa.*X.*Y) ./ ...
               (eff*beta_opt_vec.*(1-beta_opt_vec+zeta)*Y*N0 + (1-beta_opt_vec)*N0 + epsilon);
CP_Sim_OPSA = mean(0.5*log2(1+Gamma_AB_opt) > Rt, 2);

% Ground relay benchmark (Rayleigh)
X_g = d_ar^(-alpha_N) * exprnd(1, [n, N]);
Y_g = d_rb^(-alpha_N) * exprnd(1, [n, N]);
epsilon_g = zeta*N0^2 ./ (Pa.*X_g + Pb.*Y_g);
Gamma_AB_g = (eff*beta*(1-beta)*Pa.*X_g.*Y_g) ./ ...
             (eff*beta*(1-beta+zeta)*Y_g*N0 + (1-beta)*N0 + epsilon_g);
CP_Sim_GR = mean(0.5*log2(1+Gamma_AB_g) > Rt, 2);

% Analytical closed-form (series)
D_ord = 15; R_ord = 15;
CP_Theory = zeros(1, n);
for i = 1:n
    Pa_i = PVec(i)*lambda;
    A = (1-beta+zeta)*N0*delta_t / (1-beta) / Pa_i / L_au;
    B = N0*delta_t / eff / beta / Pa_i / L_au / L_ub;
    Func = 0;
    for d = 0:D_ord
        for u = 0:d
            for s = 0:u
                for r = 0:R_ord
                    Func = Func + 2 * ...
                        (gamma(d+D_ord)*D_ord^(1-2*d)*K_au^d*(1+K_au)^u * ...
                         B^(u-s)*A^s*R_ord^(1-2*r)*nchoosek(u,s)*gamma(R_ord+r)) / ...
                        (gamma(d+1)*gamma(D_ord-d+1)*gamma(u-s+1)*gamma(s+1) * ...
                         (gamma(r+1))^2*gamma(R_ord-r+1)) * ...
                        (K_ub*(1+K_ub))^r * ((1+K_au)*B/(1+K_ub))^((r-(u-s)+1)/2) * ...
                        besselk(r-(u-s)+1, 2*sqrt((1+K_au)*B*(1+K_ub)));
                end
            end
        end
    end
    CP_Theory(i) = (1+K_ub)*exp(-K_au-K_ub-(1+K_au)*A)*Func;
end

fig1 = figure('Name', 'Figure 1: Connection Probability', 'NumberTitle', 'off');
plot(P_dBW, CP_Theory,   'pr',  'MarkerSize', 8,   'LineWidth', 1, 'DisplayName', 'UR-EPSA Theory');
hold on;
plot(P_dBW, CP_Sim_EPA,  '-b',  'LineWidth', 1.5, 'DisplayName', 'UR-EPSA Simulation');
plot(P_dBW, CP_Sim_OPSA, '-og', 'LineWidth', 1.5, 'DisplayName', 'UR-OPSA Simulation');
plot(P_dBW, CP_Sim_GR,   '-k',  'LineWidth', 1.5, 'DisplayName', 'GR-EPSA Simulation');
legend('show', 'Location', 'best', 'FontSize', 10, 'FontName', 'Times New Roman');
xlabel('Network Transmit Power P [dBW]', 'FontSize', 14, 'FontName', 'Times New Roman');
ylabel('Connection Probability',          'FontSize', 14, 'FontName', 'Times New Roman');
grid on; box on;

%% ======================================================================
%% FIGURE 2: Intercept Probability vs Transmit Power
%% ======================================================================
fprintf('=== Figure 2: Intercept Probability ===\n');
P_dBW2 = linspace(0, 10, n);
PVec2  = db2pow(P_dBW2);
beta2  = 0.5;

lambda_cj = 0.7;    % With cooperative jamming
Pa_cj  = repmat(PVec2', 1, N) * lambda_cj;
Pb_cj  = repmat(PVec2', 1, N) * (1 - lambda_cj);
epsilon_cj = zeta*N0^2 ./ (Pa_cj.*X + Pb_cj.*Y);
Gamma_e1_cj = (Pa_cj/N0).*V./((Pb_cj/N0).*W + 1);
Gamma_e2_cj = eff*beta2*(1-beta2)*Pa_cj.*X.*Z ./ ...
    ((eff*beta2*(1-beta2)*Pb_cj.*Y.*Z + Z*eff*beta2*(1-beta2+zeta)*N0 + (1-beta2)*N0 + epsilon_cj));
IP_Sim_CJ = mean(max(Gamma_e1_cj, Gamma_e2_cj) > delta_e, 2);

lambda_no = 1;      % Without CJ
Pa_no  = repmat(PVec2', 1, N) * lambda_no;
Pb_no  = repmat(PVec2', 1, N) * 0;
epsilon_no = zeta*N0^2 ./ (Pa_no.*X + Pb_no.*Y);
Gamma_e1_no = (Pa_no/N0).*V ./ (1);
Gamma_e2_no = eff*beta2*(1-beta2)*Pa_no.*X.*Z ./ ...
    ((Z*eff*beta2*(1-beta2+zeta)*N0 + (1-beta2)*N0 + epsilon_no));
IP_Sim_NoCJ = mean(max(Gamma_e1_no, Gamma_e2_no) > delta_e, 2);

% Ground relay benchmark
X_g2 = d_ar^(-alpha_N)*exprnd(1,[n,N]); Y_g2 = d_rb^(-alpha_N)*exprnd(1,[n,N]);
Z_g2 = d_re^(-alpha_N)*exprnd(1,[n,N]);
epsilon_g2 = zeta*N0^2 ./ (Pa_cj.*X_g2 + Pb_cj.*Y_g2);
Gamma_e2_gr = eff*beta2*(1-beta2)*Pa_cj.*X_g2.*Z_g2 ./ ...
    ((eff*beta2*(1-beta2)*Pb_cj.*Y_g2.*Z_g2 + Z_g2*eff*beta2*(1-beta2+zeta)*N0 + ...
      (1-beta2)*N0 + epsilon_g2));
IP_Sim_GR = mean(max(Gamma_e1_cj, Gamma_e2_gr) > delta_e, 2);

fig2 = figure('Name', 'Figure 2: Intercept Probability', 'NumberTitle', 'off');
plot(P_dBW2, IP_Sim_CJ,    '-b',  'LineWidth', 1.5, 'DisplayName', 'UR - with CJ');
hold on;
plot(P_dBW2, IP_Sim_NoCJ,  '--b', 'LineWidth', 1.5, 'DisplayName', 'UR - without CJ');
plot(P_dBW2, IP_Sim_GR,    '--r', 'LineWidth', 1.5, 'DisplayName', 'GR - without CJ');
legend('show', 'Location', 'best', 'FontSize', 10, 'FontName', 'Times New Roman');
xlabel('Network Transmit Power P [dBW]', 'FontSize', 14, 'FontName', 'Times New Roman');
ylabel('Intercept Probability',           'FontSize', 14, 'FontName', 'Times New Roman');
grid on; box on;

%% ======================================================================
%% FIGURE 3: Average Secrecy Rate vs lambda and beta (3D surface)
%% ======================================================================
fprintf('=== Figure 3: ESR Surface vs lambda and beta ===\n');
betaVec2  = linspace(0, 1, n);
lambdaVec = linspace(0, 1, n);
Rs_Surface = zeros(n, n);
P_fixed = db2pow(20);

for i = 1:n
    lambda_i = lambdaVec(i);
    beta_i   = betaVec2;  % row sweep
    Pa_i = lambda_i * P_fixed;
    Pb_i = (1 - lambda_i) * P_fixed;
    X1 = X(1, :);  Y1 = Y(1, :);  Z1 = Z(1, :);
    V1 = V(1, :);  W1 = W(1, :);
    epsilon_i = zeta*N0^2 ./ (Pa_i*X1 + Pb_i*Y1);

    for j = 1:n
        b_ij = betaVec2(j);
        G_AB = (eff*b_ij*(1-b_ij)*Pa_i .* X1 .* Y1) ./ ...
               (eff*b_ij*(1-b_ij+zeta)*Y1*N0 + (1-b_ij)*N0 + epsilon_i);
        G_e1 = (Pa_i/N0).*V1./((Pb_i/N0).*W1 + 1);
        G_e2 = eff*b_ij*(1-b_ij)*Pa_i.*X1.*Z1 ./ ...
               ((eff*b_ij*(1-b_ij)*Pb_i.*Y1.*Z1 + Z1*eff*b_ij*(1-b_ij+zeta)*N0 + ...
                 (1-b_ij)*N0 + epsilon_i));
        G_E = max(G_e1, G_e2);
        Rs_Surface(i, j) = mean(max(0.5*log2(1+G_AB) - 0.5*log2(1+G_E), 0));
    end
end

fig3 = figure('Name', 'Figure 3: ESR Surface', 'NumberTitle', 'off');
[LV, BV] = meshgrid(lambdaVec, betaVec2);
surf(LV, BV, Rs_Surface', 'EdgeAlpha', 0.3);
colorbar;
xlabel('Power Allocation Factor (\lambda)', 'FontSize', 14, 'FontName', 'Times New Roman');
ylabel('Power Splitting Ratio (\beta)',     'FontSize', 14, 'FontName', 'Times New Roman');
zlabel('Avg Secrecy Rate [bits/s/Hz]',      'FontSize', 14, 'FontName', 'Times New Roman');
title('Average Secrecy Rate Surface', 'FontName', 'Times New Roman', 'FontSize', 13);
grid on; view(45, 30);

%% ========== Save All Figures ============================================
results_dir = fullfile(fileparts(mfilename('fullpath')), '..', '..', 'results', 'figures');
if ~exist(results_dir, 'dir'); mkdir(results_dir); end
savefig(fig1, fullfile(results_dir, 'Fig1_CP_vs_Power.fig'));
print(fig1,   fullfile(results_dir, 'Fig1_CP_vs_Power.png'), '-dpng', '-r300');
savefig(fig2, fullfile(results_dir, 'Fig2_IP_vs_Power.fig'));
print(fig2,   fullfile(results_dir, 'Fig2_IP_vs_Power.png'), '-dpng', '-r300');
savefig(fig3, fullfile(results_dir, 'Fig3_ESR_Surface.fig'));
print(fig3,   fullfile(results_dir, 'Fig3_ESR_Surface.png'), '-dpng', '-r300');

data_dir = fullfile(fileparts(mfilename('fullpath')), '..', '..', 'data', 'processed');
if ~exist(data_dir, 'dir'); mkdir(data_dir); end
save(fullfile(data_dir, 'Comprehensive_Results.mat'), ...
     'P_dBW', 'CP_Sim_EPA', 'CP_Sim_OPSA', 'CP_Sim_GR', 'CP_Theory', ...
     'P_dBW2', 'IP_Sim_CJ', 'IP_Sim_NoCJ', 'IP_Sim_GR', ...
     'lambdaVec', 'betaVec2', 'Rs_Surface', 'N');

fprintf('\nAll figures saved to results/figures/\n');
fprintf('All data saved to data/processed/Comprehensive_Results.mat\n');

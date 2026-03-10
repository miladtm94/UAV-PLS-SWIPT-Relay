%% INTERCEPT PROBABILITY vs. NETWORK TRANSMIT POWER
%
%  Computes the probability that the eavesdropper can decode the secret
%  message (i.e., the eavesdropper's SNR exceeds a threshold gamma0).
%
%  IP = Pr{max(Gamma_e1, Gamma_e2) > gamma0}
%     = 1 - Pr{Gamma_e1 < gamma0} * Pr{Gamma_e2 < gamma0}
%
%  Compares simulation with 2D numerical integration analytical result.
%
%  Reference: Fig. X in the published paper.
%
% -------------------------------------------------------------------------

clc; clear; close all;

addpath(fullfile(fileparts(mfilename('fullpath')), '..', '..', 'src', 'analytics'));

%% ========== Network Topology ============================================
D = 5; H = 1;
t_A = [0, 0, 0]; t_U = [D/5, 0, H]; t_B = [D, 0, 0]; t_E = [4*D/5, 1, 0];
d_au = norm(t_A - t_U); d_ub = norm(t_U - t_B);
d_ae = norm(t_A - t_E); d_be = norm(t_B - t_E); d_ue = norm(t_U - t_E);

%% ========== Channel Parameters ==========================================
alpha_L = 2; alpha_N = 3.5; N0 = 1e-2; eff = 0.7; zeta = 2; beta = 0.7;
a = 0.28; b = 9.61; km = 1; kM = 7;

theta_ub = asin(H/d_ub); P_los_ub = 1/(1+a*exp(-b*(theta_ub-a)));
L_ub = d_ub^(-((alpha_L-alpha_N)*P_los_ub+alpha_N));
K_ub = db2pow(km+(kM-km)*2*theta_ub/pi);

theta_au = asin(H/d_au); P_los_au = 1/(1+a*exp(-b*(theta_au-a)));
L_au = d_au^(-((alpha_L-alpha_N)*P_los_au+alpha_N));
K_au = db2pow(km+(kM-km)*2*theta_au/pi);

theta_ue = asin(H/d_ue); P_los_ue = 1/(1+a*exp(-b*(theta_ue-a)));
L_ue = d_ue^(-((alpha_L-alpha_N)*P_los_ue+alpha_N));
K_ue = db2pow(km+(kM-km)*2*theta_ue/pi);

L_ae = d_ae^(-alpha_N); L_be = d_be^(-alpha_N);

%% ========== Simulation Setup ============================================
n = 5; N = 1e5; gamma0 = 1;
P_dBW = linspace(0, 25, n);
P_lin = db2pow(P_dBW);
PVec  = repmat(P_lin', 1, N);
Pa = PVec/2; Pb = PVec/2;

X = ncx2rnd(2,sqrt(2*K_au),[n,N]); X=X./mean(X); X=X*L_au;
Y = ncx2rnd(2,sqrt(2*K_ub),[n,N]); Y=Y./mean(Y); Y=Y*L_ub;
Z = ncx2rnd(2,sqrt(2*K_ue),[n,N]); Z=Z./mean(Z); Z=Z*L_ue;
V = L_ae*exprnd(1,[n,N]); W = L_be*exprnd(1,[n,N]);

%% ========== Monte Carlo Simulation =====================================
epsilon  = zeta*N0^2 ./ (Pa.*X + Pb.*Y);
Gamma_e1 = (Pa/N0).*V ./ ((Pb/N0).*W + 1);
Gamma_e2 = eff*beta*(1-beta)*(Pa/N0).*X.*Z ./ ...
           (eff*beta*(1-beta)*(Pb/N0).*Y.*Z + Z*eff*beta*(1-beta+zeta) + 1-beta);
Gamma_E  = max(Gamma_e1, Gamma_e2);
IP_Sim   = mean(Gamma_E > gamma0, 2);

%% ========== Analytical (2D Numerical Integration) ======================
P_gammaE1 = zeros(1, n);
P_gammaE2 = zeros(1, n);
IP_Analytic = zeros(1, n);

y_grid = [linspace(1e-4, 1, 150), linspace(1, 10, 100)];
z_grid = [linspace(1e-4, 1, 150), linspace(1, 10, 100)];
[YY, ZZ] = meshgrid(y_grid, z_grid);
fZ = (K_ue+1)*exp(-K_ue)*exp(-(K_ue+1)*ZZ).*besseli(0,2*sqrt(K_ue*(K_ue+1)*ZZ));
fY = (K_ub+1)*exp(-K_ub)*exp(-(K_ub+1)*YY).*besseli(0,2*sqrt(K_ub*(K_ub+1)*YY));

for i = 1:n
    Pa_i = P_lin(i)/2; Pb_i = P_lin(i)/2;
    a1 = Pb_i*L_ub / Pa_i / L_au;
    a2 = N0 / eff / beta / Pa_i / L_au / L_ue;
    a3 = (1-beta+zeta)*N0 / (1-beta) / Pa_i / L_au;
    c  = 2*(1+K_au)*gamma0*(a2./ZZ + a3);
    bb = 2*(1+K_au)*gamma0*a1;

    func = marcumq(sqrt(2*K_au), sqrt(bb*YY + c)) .* fY .* fZ;
    P_gammaE2(i) = trapz(z_grid, trapz(y_grid, func, 2));
    P_gammaE1(i) = 1 - Pa_i*L_ae*exp(-N0*gamma0/Pa_i/L_ae) / (Pb_i*L_be*gamma0 + Pa_i*L_ae);
    IP_Analytic(i) = 1 - P_gammaE1(i) * P_gammaE2(i);
end

%% ========== Plot ========================================================
figure('Name', 'Intercept Probability vs Transmit Power');
plot(P_dBW, IP_Sim,     '-b',  'LineWidth', 1.5, 'DisplayName', 'Simulation');
hold on;
plot(P_dBW, IP_Analytic,'pr',  'MarkerSize', 8,  'DisplayName', 'Analytical');
xlabel('Network Transmit Power P [dBW]', 'FontSize', 14, 'FontName', 'Times New Roman');
ylabel('Intercept Probability',           'FontSize', 14, 'FontName', 'Times New Roman');
legend('show', 'Location', 'best', 'FontSize', 11);
grid on; box on;

%% ========== Save ========================================================
results_dir = fullfile(fileparts(mfilename('fullpath')), '..', '..', 'results', 'figures');
if ~exist(results_dir, 'dir'); mkdir(results_dir); end
savefig(fullfile(results_dir, 'IP_vs_Power.fig'));
print( fullfile(results_dir, 'IP_vs_Power.png'), '-dpng', '-r300');

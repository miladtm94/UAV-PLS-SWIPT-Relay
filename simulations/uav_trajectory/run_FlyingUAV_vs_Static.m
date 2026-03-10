%% FLYING UAV RELAY vs. STATIC GROUND RELAY
%
%  Compares Average Secrecy Rate (ASR), Intercept Probability (IP), and
%  Connection Probability (CP) for:
%    (a) Flying UAV-based relay with Equal Power Allocation  (EPA)
%    (b) Flying UAV-based relay with Optimal Power Allocation (OPA)
%    (c) Static ground relay with Equal Power Allocation
%
%  The UAV moves linearly from Alice's location to Bob's location at
%  constant velocity v. At each time step, the UAV position is updated and
%  instantaneous metrics are averaged over the full trajectory.
%
%  Reference: Fig. X in the published paper.
%
% -------------------------------------------------------------------------

clc; clear; close all;

addpath(fullfile(fileparts(mfilename('fullpath')), '..', '..', 'src', 'analytics'));
addpath(fullfile(fileparts(mfilename('fullpath')), '..', '..', 'src', 'channel'));

%% ========== Scenario Parameters =========================================
D    = 1000;   % Alice-to-Bob distance (m)
H    = 100;    % UAV flight altitude (m)
v    = 50;     % UAV speed (m/s)
n    = D / v;  % Number of time steps (= 20 steps)

% Eavesdropper location
t_A  = [0,       0, 0];
t_B  = [D,       0, 0];
t_E  = [4*D/5,   0, 0];

%% ========== Channel & SWIPT Parameters ==================================
alpha_L = 2; alpha_N = 4;
N0  = 1e-15; eff = 0.7; zeta = 2; beta = 0.7;

% A2G LoS probability parameters (urban)
a = 0.28; b = 9.61;
km = 5; kM = 15;   % Rician K-factor range [dB]

N = 1e5;   % Monte Carlo samples

% Fixed transmit power
P_dBW = 20;
P     = db2pow(P_dBW);

%% ========== UAV Trajectory: EPA =========================================
fprintf('Simulating flying UAV relay (EPA)...\n');
Rs_UAV_EPA = zeros(1, n + 1);
IP_UAV_EPA = zeros(1, n + 1);
CP_UAV_EPA = zeros(1, n + 1);
Rt = 2; Re = 1;

lambda_epa = 0.5;
Pa_epa = lambda_epa * P;
Pb_epa = (1 - lambda_epa) * P;

for i = 1:n + 1
    t_U = [v*(i-1), 0, H];
    [L_au, K_au] = compute_a2g_stats(t_A, t_U, alpha_L, alpha_N, a, b, km, kM);
    [L_ub, K_ub] = compute_a2g_stats(t_U, t_B, alpha_L, alpha_N, a, b, km, kM);
    [L_ue, K_ue] = compute_a2g_stats(t_U, t_E, alpha_L, alpha_N, a, b, km, kM);
    d_ae = norm(t_A - t_E);  L_ae = d_ae^(-alpha_N);
    d_be = norm(t_B - t_E);  L_be = d_be^(-alpha_N);

    [X, ~] = generate_rician(K_au, L_au, 1, N);
    [Y, ~] = generate_rician(K_ub, L_ub, 1, N);
    [Z, ~] = generate_rician(K_ue, L_ue, 1, N);
    V = L_ae * exprnd(1, [1, N]);
    W = L_be * exprnd(1, [1, N]);

    epsilon  = zeta * N0^2 ./ (Pa_epa * X + Pb_epa * Y);
    Gamma_AB = (eff*beta*(1-beta)*Pa_epa .* X .* Y) ./ ...
               (eff*beta*(1-beta+zeta)*Y*N0 + (1-beta)*N0 + epsilon);
    Gamma_e1 = (Pa_epa/N0) .* V ./ ((Pb_epa/N0) .* W + 1);
    Gamma_e2 = eff*beta*(1-beta)*Pa_epa .* X .* Z ./ ...
               ((eff*beta*(1-beta)*Pb_epa.*Y.*Z + Z*eff*beta*(1-beta+zeta)*N0 + (1-beta)*N0 + epsilon));
    Gamma_E  = max(Gamma_e1, Gamma_e2);
    Cs = max(0.5*log2(1+Gamma_AB) - 0.5*log2(1+Gamma_E), 0);
    Rs_UAV_EPA(i) = mean(Cs);
    IP_UAV_EPA(i) = mean(0.5*log2(1+Gamma_E) > Re);
    CP_UAV_EPA(i) = mean(0.5*log2(1+Gamma_AB) > Rt);
end

%% ========== UAV Trajectory: OPA =========================================
fprintf('Simulating flying UAV relay (OPA)...\n');
Rs_UAV_OPA = zeros(1, n + 1);
IP_UAV_OPA = zeros(1, n + 1);
CP_UAV_OPA = zeros(1, n + 1);

for i = 1:n + 1
    t_U = [v*(i-1), 0, H];
    [L_au, K_au] = compute_a2g_stats(t_A, t_U, alpha_L, alpha_N, a, b, km, kM);
    [L_ub, K_ub] = compute_a2g_stats(t_U, t_B, alpha_L, alpha_N, a, b, km, kM);
    [L_ue, K_ue] = compute_a2g_stats(t_U, t_E, alpha_L, alpha_N, a, b, km, kM);
    d_ae = norm(t_A - t_E);  L_ae = d_ae^(-alpha_N);
    d_be = norm(t_B - t_E);  L_be = d_be^(-alpha_N);

    [X, ~] = generate_rician(K_au, L_au, 1, N);
    [Y, ~] = generate_rician(K_ub, L_ub, 1, N);
    [Z, ~] = generate_rician(K_ue, L_ue, 1, N);
    V = L_ae * exprnd(1, [1, N]);
    W = L_be * exprnd(1, [1, N]);

    % Optimal power allocation (closed-form)
    lambda_opt = (-W.*Y + sqrt(W.*Y.*(V.*Y + W.*X))) ./ ((X - Y).*W + V.*Y);
    lambda_opt = min(max(lambda_opt, 0), 1);   % Clip to [0,1]
    Pa_opt = lambda_opt * P;
    Pb_opt = (1 - lambda_opt) * P;

    epsilon  = zeta * N0^2 ./ (Pa_opt .* X + Pb_opt .* Y);
    Gamma_AB = (eff*beta*(1-beta)*Pa_opt .* X .* Y) ./ ...
               (eff*beta*(1-beta+zeta)*Y*N0 + (1-beta)*N0 + epsilon);
    Gamma_e1 = (Pa_opt/N0) .* V ./ ((Pb_opt/N0) .* W + 1);
    Gamma_e2 = eff*beta*(1-beta)*Pa_opt .* X .* Z ./ ...
               ((eff*beta*(1-beta)*Pb_opt.*Y.*Z + Z*eff*beta*(1-beta+zeta)*N0 + (1-beta)*N0 + epsilon));
    Gamma_E  = max(Gamma_e1, Gamma_e2);
    Cs = max(0.5*log2(1+Gamma_AB) - 0.5*log2(1+Gamma_E), 0);
    Rs_UAV_OPA(i) = mean(Cs);
    IP_UAV_OPA(i) = mean(0.5*log2(1+Gamma_E) > Re);
    CP_UAV_OPA(i) = mean(0.5*log2(1+Gamma_AB) > Rt);
end

%% ========== Static Ground Relay (EPA) ===================================
fprintf('Simulating static ground relay (EPA)...\n');
Rs_Static = zeros(1, n + 1);
IP_Static = zeros(1, n + 1);
CP_Static = zeros(1, n + 1);

t_R_static = [D/2, 0, 0];   % Midpoint static relay
d_ar = norm(t_A - t_R_static);  L_ar = d_ar^(-alpha_N);
d_rb = norm(t_R_static - t_B);  L_rb = d_rb^(-alpha_N);
d_re = norm(t_R_static - t_E);  L_re = d_re^(-alpha_N);
d_ae = norm(t_A - t_E);  L_ae = d_ae^(-alpha_N);
d_be = norm(t_B - t_E);  L_be = d_be^(-alpha_N);

Pa_s = P/2;  Pb_s = P/2;
X = L_ar * exprnd(1, [1, N]);   Y = L_rb * exprnd(1, [1, N]);
Z = L_re * exprnd(1, [1, N]);   V = L_ae * exprnd(1, [1, N]);
W = L_be * exprnd(1, [1, N]);

epsilon  = zeta * N0^2 ./ (Pa_s * X + Pb_s * Y);
Gamma_AB = (eff*beta*(1-beta)*Pa_s .* X .* Y) ./ ...
           (eff*beta*(1-beta+zeta)*Y*N0 + (1-beta)*N0 + epsilon);
Gamma_e1 = (Pa_s/N0) .* V ./ ((Pb_s/N0) .* W + 1);
Gamma_e2 = eff*beta*(1-beta)*Pa_s .* X .* Z ./ ...
           ((eff*beta*(1-beta)*Pb_s.*Y.*Z + Z*eff*beta*(1-beta+zeta)*N0 + (1-beta)*N0 + epsilon));
Gamma_E  = max(Gamma_e1, Gamma_e2);
Cs = max(0.5*log2(1+Gamma_AB) - 0.5*log2(1+Gamma_E), 0);
Rs_Static_val = mean(Cs);
IP_Static_val = mean(0.5*log2(1+Gamma_E) > Re);
CP_Static_val = mean(0.5*log2(1+Gamma_AB) > Rt);
Rs_Static(:) = Rs_Static_val;
IP_Static(:) = IP_Static_val;
CP_Static(:) = CP_Static_val;

%% ========== Plot ========================================================
x_pos = v * (0:n);

figure('Name', 'Flying UAV Relay vs Static Relay', 'Position', [100 100 1200 350]);

subplot(1, 3, 1);
plot(x_pos, Rs_UAV_OPA, '-sg', 'LineWidth', 1.5, 'DisplayName', 'UAV Relay - OPA');
hold on;
plot(x_pos, Rs_UAV_EPA, '-or', 'LineWidth', 1.5, 'DisplayName', 'UAV Relay - EPA');
plot(x_pos, Rs_Static,  '--b', 'LineWidth', 1.5, 'DisplayName', 'Static Relay - EPA');
legend('show', 'Location', 'best', 'FontSize', 9);
xlabel('UAV x-position [m]', 'FontName', 'Times New Roman', 'FontSize', 12);
ylabel('Average Secrecy Rate [bits/s/Hz]', 'FontName', 'Times New Roman', 'FontSize', 12);
grid on;

subplot(1, 3, 2);
plot(x_pos, IP_UAV_OPA, '-sg', 'LineWidth', 1.5, 'DisplayName', 'UAV Relay - OPA');
hold on;
plot(x_pos, IP_UAV_EPA, '-or', 'LineWidth', 1.5, 'DisplayName', 'UAV Relay - EPA');
plot(x_pos, IP_Static,  '--b', 'LineWidth', 1.5, 'DisplayName', 'Static Relay - EPA');
xlabel('UAV x-position [m]', 'FontName', 'Times New Roman', 'FontSize', 12);
ylabel('Intercept Probability', 'FontName', 'Times New Roman', 'FontSize', 12);
legend('show', 'Location', 'best', 'FontSize', 9);
grid on;

subplot(1, 3, 3);
plot(x_pos, CP_UAV_OPA, '-sg', 'LineWidth', 1.5, 'DisplayName', 'UAV Relay - OPA');
hold on;
plot(x_pos, CP_UAV_EPA, '-or', 'LineWidth', 1.5, 'DisplayName', 'UAV Relay - EPA');
plot(x_pos, CP_Static,  '--b', 'LineWidth', 1.5, 'DisplayName', 'Static Relay - EPA');
xlabel('UAV x-position [m]', 'FontName', 'Times New Roman', 'FontSize', 12);
ylabel('Connection Probability', 'FontName', 'Times New Roman', 'FontSize', 12);
legend('show', 'Location', 'best', 'FontSize', 9);
grid on;

sgtitle('Flying UAV Relay vs Static Ground Relay', 'FontName', 'Times New Roman', 'FontSize', 14);

%% ========== Save ========================================================
results_dir = fullfile(fileparts(mfilename('fullpath')), '..', '..', 'results', 'figures');
if ~exist(results_dir, 'dir'); mkdir(results_dir); end
savefig(fullfile(results_dir, 'FlyingUAV_vs_StaticRelay.fig'));
print( fullfile(results_dir, 'FlyingUAV_vs_StaticRelay.png'), '-dpng', '-r300');


%% ========== Local Helper Functions ======================================
function [L, K] = compute_a2g_stats(t_node1, t_node2, alpha_L, alpha_N, a, b, km, kM)
    d_horiz = norm(t_node1(1:2) - t_node2(1:2));
    H       = abs(t_node2(3) - t_node1(3));
    d_3D    = norm(t_node1 - t_node2);
    if d_3D < 1e-9; L = 1; K = db2pow(kM); return; end
    theta   = asin(H / d_3D);
    P_los   = 1 / (1 + a * exp(-b * (theta - a)));
    alpha   = (alpha_L - alpha_N) * P_los + alpha_N;
    L       = d_3D^(-alpha);
    K       = db2pow(km + (kM - km) * 2 * theta / pi);
end

function [X, lam] = generate_rician(K, L, nr, nc)
    lam   = sqrt(2 * K);
    X_raw = ncx2rnd(2, lam, [nr, nc]);
    X     = (X_raw ./ mean(X_raw(:))) * L;
end

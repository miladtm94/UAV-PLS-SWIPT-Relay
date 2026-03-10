%% QUICK START EXAMPLE — UAV-PLS-SWIPT-Relay
%
%  This script demonstrates the core simulation pipeline of the paper in a
%  compact form. It computes the Average Secrecy Rate for a UAV relay network
%  with SWIPT and cooperative jamming, and plots ESR vs transmit power.
%
%  Run time: ~30 seconds (N = 1e4 samples).
%
%  Prerequisites: MATLAB R2019b or later, Statistics & Machine Learning Toolbox.
%
% -------------------------------------------------------------------------

clc; clear; close all;

% Add source paths
repo_root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(repo_root, 'src', 'analytics'));

fprintf('=== UAV-PLS-SWIPT-Relay: Quick Start Example ===\n\n');

%% ========== System Parameters ===========================================
D = 10;   H = 1.5;             % Topology (normalized)
t_A = [0, 0, 0]; t_U = [D/5, 0, H]; t_B = [D, 0, 0]; t_E = [4*D/5, 1, 0];
d_au = norm(t_A-t_U); d_ub = norm(t_U-t_B);
d_ae = norm(t_A-t_E); d_be = norm(t_B-t_E); d_ue = norm(t_U-t_E);

alpha_L = 2; alpha_N = 3.5;
N0 = 1e-2; eff = 0.7; zeta = 2;
a = 0.28; b = 9.61; km = 1; kM = 10;
beta   = 0.5;    % Power splitting ratio
lambda = 0.7;    % Power allocation factor

% A2G channel statistics
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

%% ========== Monte Carlo ESR =============================================
n = 20; N = 1e4;
P_dBW = linspace(5, 30, n);
P_lin = db2pow(P_dBW);
ESR_sim = zeros(1, n);

X = ncx2rnd(2,sqrt(2*K_au),[n,N]); X=X./mean(X); X=X*L_au;
Y = ncx2rnd(2,sqrt(2*K_ub),[n,N]); Y=Y./mean(Y); Y=Y*L_ub;
Z = ncx2rnd(2,sqrt(2*K_ue),[n,N]); Z=Z./mean(Z); Z=Z*L_ue;
V = L_ae*exprnd(1,[n,N]); W = L_be*exprnd(1,[n,N]);

for i = 1:n
    Pa = lambda * P_lin(i); Pb = (1-lambda) * P_lin(i);
    epsilon = zeta*N0^2 ./ (Pa.*X(i,:) + Pb.*Y(i,:));
    Gamma_AB = (eff*beta*(1-beta)*Pa.*X(i,:).*Y(i,:)) ./ ...
               (eff*beta*(1-beta+zeta)*Y(i,:)*N0 + (1-beta)*N0 + epsilon);
    Gamma_e1 = (Pa/N0).*V(i,:) ./ ((Pb/N0).*W(i,:) + 1);
    Gamma_e2 = eff*beta*(1-beta)*Pa.*X(i,:).*Z(i,:) ./ ...
               ((eff*beta*(1-beta)*Pb.*Y(i,:).*Z(i,:) + Z(i,:)*eff*beta*(1-beta+zeta)*N0 ...
                 + (1-beta)*N0 + epsilon));
    Gamma_E = max(Gamma_e1, Gamma_e2);
    ESR_sim(i) = mean(max(0.5*log2(1+Gamma_AB) - 0.5*log2(1+Gamma_E), 0));
end

%% ========== Plot ========================================================
figure('Name', 'Quick Start: ESR vs Transmit Power');
plot(P_dBW, ESR_sim, '-or', 'LineWidth', 2, 'MarkerFaceColor', 'r');
xlabel('Network Transmit Power P [dBW]', 'FontSize', 14, 'FontName', 'Times New Roman');
ylabel('Average Secrecy Rate [bits/s/Hz]', 'FontSize', 14, 'FontName', 'Times New Roman');
title(sprintf('ESR vs Power (\\beta=%.1f, \\lambda=%.1f)', beta, lambda), ...
      'FontSize', 13, 'FontName', 'Times New Roman');
grid on; box on;

fprintf('ESR at P=20 dBW: %.4f bits/s/Hz\n', interp1(P_dBW, ESR_sim, 20));
fprintf('\nTo reproduce all paper figures, see:\n');
fprintf('  simulations/comprehensive/run_Comprehensive_Simulation.m\n');
fprintf('  simulations/connection_probability/run_CP_vs_Power.m\n');
fprintf('  simulations/secrecy_rate/run_ESR_vs_SNR.m\n');
fprintf('  simulations/uav_trajectory/run_FlyingUAV_vs_Static.m\n');

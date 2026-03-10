%% SECRECY OUTAGE PROBABILITY vs. TRANSMIT POWER
%
%  Computes the Secrecy Outage Probability (SOP) for a UAV relay network
%  with multiple randomly located eavesdroppers following a Poisson Point
%  Process (PPP) within a disk of radius Re centered at Bob.
%
%  SOP = Probability that the instantaneous secrecy rate < target Rs.
%
%  Compares: Monte Carlo simulation vs. analytical expression derived
%  in the paper (which leverages the PGFL of the PPP).
%
%  Reference: Fig. X in the published paper.
%
% -------------------------------------------------------------------------

clc; clear; close all;

addpath(fullfile(fileparts(mfilename('fullpath')), '..', '..', 'src', 'analytics'));

%% ========== Network Topology ============================================
D  = 5;   H = 50;
d1 = D/3; d2 = D/2;   % UAV horizontal offsets
d_au = sqrt(d1^2 + H^2);
d_ub = sqrt(d2^2 + H^2);

% Eavesdropper PPP model
Re      = 30;         % Radius of eavesdropper disk around Bob (m)
area    = pi * Re^2;
lambda_E = 0.001;     % PPP intensity [eavesdroppers/m^2]

%% ========== Channel Parameters ==========================================
a_los  = 2; a_nlos = 3;
eta_los_dB = 1;  eta_los  = db2pow(eta_los_dB);
eta_nlos_dB = 20; eta_nlos = db2pow(eta_nlos_dB);
N0_dB = -70;  N0 = 1e-3 * db2pow(N0_dB);
e     = 0.5;   alpha = 0.5;   R_Loss = 0.5;
delta_e = 2^(2*R_Loss) - 1;

% Non-centrality parameter (empirical model)
lambda_au = -17.4 * log10(H) + 29.6;
lambda_ue = -17.4 * log10(H) + 29.6;

L_au = 1 / (eta_los * d_au^a_los);

% LoS model constants
a_prob = 9.61; b_prob = 0.28;

%% ========== Simulation Setup ============================================
n  = 20;
Ne = 1e3;  % Outer Monte Carlo (PPP realizations)
N  = 1e3;  % Inner Monte Carlo (channel)
Pa_dBW = linspace(0, 25, n);
Pso_Simulation = zeros(1, n);
Pso_Analytic   = zeros(1, n);

%% ========== Monte Carlo Simulation =====================================
fprintf('Running Monte Carlo simulation (PPP eavesdroppers)...\n');
for i = 1:n
    Pa  = db2pow(Pa_dBW(i));
    rho = Pa / N0;
    X   = ncx2rnd(2, lambda_au, [1, N]);
    G_outer = zeros(1, Ne);

    for k = 1:Ne
        NumEs  = poissrnd(area * lambda_E);
        PhiE   = 2*pi * rand(NumEs, 1);
        rE     = Re * sqrt(rand(NumEs, 1));
        F      = 1;

        for j = 1:NumEs
            re = rE(j); Phie = PhiE(j);
            d_ue = sqrt(H^2 + d2^2 + re^2 - 2*d2*re*cos(Phie));
            P_los_ue  = 1 / (1 + a_prob * exp(-b_prob * (asin(H/d_ue) - a_prob)));
            L_los_ue  = 1 / (eta_los  * d_ue^a_los);
            L_nlos_ue = 1 / (eta_nlos * d_ue^a_nlos);
            L_ue      = P_los_ue * L_los_ue + (1 - P_los_ue) * L_nlos_ue;
            Y = ncx2rnd(2, lambda_ue, [1, N]);
            Gamma_e = (e*alpha*(1-alpha)*rho^2 * L_au^2 * L_ue .* X.^2 .* Y) ./ ...
                      (e*alpha*rho*L_au*L_ue.*X.*Y + (1-alpha)*rho*L_au.*X + 1);
            F = F * mean(Gamma_e < delta_e);
        end
        G_outer(k) = F;
    end
    Pso_Simulation(i) = 1 - mean(G_outer);
end

%% ========== Analytical Expression =======================================
fprintf('Computing analytical SOP...\n');
D_ord = 1; R_ord = 1;

for i = 1:n
    Pa  = db2pow(Pa_dBW(i));
    rho = Pa / N0;
    Ae  = delta_e / (1 - alpha) / rho / L_au;
    L   = 0;

    for d = 0:D_ord
        for u = 0:d
            for v = 0:u
                for r = 0:R_ord
                    Num = gamma(d+D_ord) * D_ord^(1-2*d) * lambda_au^d * factorial(d) ...
                          * lambda_ue^r * Ae^(u-v) * R_ord^(1-2*r) * nchoosek(u,v) ...
                          * gamma(R_ord+r) * 2^(-u-2*r-d);
                    Den = (gamma(d+1))^2 * gamma(D_ord-d+1) * factorial(u) ...
                          * (gamma(r+1))^2 * gamma(R_ord-r+1) ...
                          * exp(0.5*(lambda_au + lambda_ue + Ae));
                    P_coef = Num / Den;
                    nu1 = (r + v + 1) * 0.5;
                    nu2 = r - v + 1;

                    d_ue_fn   = @(xe, ye) sqrt(H^2 + d2^2 + (xe.^2 + ye.^2) - 2*xe*d2);
                    L_los_fn  = @(x) 1 ./ (eta_los  .* x.^a_los);
                    P_los_fn  = @(x) 1 ./ (1 + a_prob * exp(-b_prob*(asin(H./x)-a_prob)));
                    L_nlos_fn = @(x) 1 ./ (eta_nlos .* x.^a_nlos);
                    L_ue_fn   = @(x) P_los_fn(x).*L_los_fn(x) + (1-P_los_fn(x)).*L_nlos_fn(x);
                    Be_fn     = @(xe, ye) delta_e / alpha / e / rho / L_au ./ L_ue_fn(d_ue_fn(xe, ye));
                    inner_fn  = @(xe, ye) (Be_fn(xe, ye)).^nu1 .* besselk(nu2, sqrt(Be_fn(xe, ye)));
                    polar_fn  = @(rr, th) inner_fn(rr.*cos(th), rr.*sin(th)) .* rr;
                    L = L + P_coef * integral2(polar_fn, 0, Re, 0, 2*pi);
                end
            end
        end
    end
    Pso_Analytic(i) = 1 - exp(-lambda_E * L);
end

%% ========== Plot ========================================================
figure('Name', 'Secrecy Outage Probability vs Transmit Power');
semilogy(Pa_dBW, Pso_Simulation, '-sr', 'LineWidth', 1.5, 'DisplayName', 'Simulation');
hold on;
semilogy(Pa_dBW, Pso_Analytic,  '-ob', 'LineWidth', 1.5, 'DisplayName', 'Analytical');
xlabel('Source Transmit Power P_a [dBW]', 'FontSize', 14, 'FontName', 'Times New Roman');
ylabel('Secrecy Outage Probability',       'FontSize', 14, 'FontName', 'Times New Roman');
legend('show', 'Location', 'best', 'FontSize', 11);
grid on; box on;

%% ========== Save ========================================================
results_dir = fullfile(fileparts(mfilename('fullpath')), '..', '..', 'results', 'figures');
if ~exist(results_dir, 'dir'); mkdir(results_dir); end
savefig(fullfile(results_dir, 'SOP_vs_Power.fig'));
print( fullfile(results_dir, 'SOP_vs_Power.png'), '-dpng', '-r300');

MSE = mean(abs(Pso_Simulation - Pso_Analytic) ./ max(Pso_Simulation, 1e-10));
fprintf('Mean Absolute Relative Error: %.3f%%\n', MSE * 100);

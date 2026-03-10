%% CONNECTION PROBABILITY: OPTIMAL vs EQUAL POWER ALLOCATION
%
%  Compares the Connection Probability (CP) under:
%    - Optimal Power Allocation (OPA):  lambda = 1  (all power to Alice->UAV)
%    - Equal Power Allocation  (EPA):   lambda = 0.5
%
%  Both simulation and closed-form approximation results are plotted.
%  This figure demonstrates the advantage of OPA over EPA for the
%  connection probability metric.
%
%  Reference: Fig. X in the published paper.
%
% -------------------------------------------------------------------------

clc; clear; close all;

addpath(fullfile(fileparts(mfilename('fullpath')), '..', '..', 'src', 'analytics'));

%% ========== Network Topology ============================================
D = 5;  H = 1;
t_A = [0, 0, 0]; t_U = [D/5, 0, H]; t_B = [D, 0, 0]; t_E = [4*D/5, 1, 0];
d_au = norm(t_A - t_U);  d_ub = norm(t_U - t_B);
d_ae = norm(t_A - t_E);  d_be = norm(t_B - t_E);  d_ue = norm(t_U - t_E);

%% ========== Channel Parameters ==========================================
alpha_L = 2; alpha_N = 4; a = 0.28; b = 9.61; km = 1; kM = 7;
N0 = 1e-2; eff = 0.7; zeta = 2; beta = 0.7;
Rt = 1; delta_t = 2^(2*Rt) - 1;

% UAV link statistics
theta_ub   = asin(H / d_ub);
P_los_ub   = 1 / (1 + a * exp(-b * (theta_ub - a)));
alpha_ub   = (alpha_L - alpha_N) * P_los_ub + alpha_N;
L_ub       = d_ub^(-alpha_ub);
K_ub       = db2pow(km + (kM - km) * 2 * theta_ub / pi);

theta_au   = asin(H / d_au);
P_los_au   = 1 / (1 + a * exp(-b * (theta_au - a)));
alpha_au   = (alpha_L - alpha_N) * P_los_au + alpha_N;
L_au       = d_au^(-alpha_au);
K_au       = db2pow(km + (kM - km) * 2 * theta_au / pi);

theta_ue   = asin(H / d_ue);
P_los_ue   = 1 / (1 + a * exp(-b * (theta_ue - a)));
alpha_ue   = (alpha_L - alpha_N) * P_los_ue + alpha_N;
L_ue       = d_ue^(-alpha_ue);
K_ue       = db2pow(km + (kM - km) * 2 * theta_ue / pi);

L_ae = d_ae^(-alpha_N);  L_be = d_be^(-alpha_N);

n = 20; N = 1e5;
P_dBW = linspace(0, 25, n);
PVec  = db2pow(P_dBW);
PMat  = repmat(PVec', 1, N);

%% ========== Channel Simulation ==========================================
X = ncx2rnd(2, sqrt(2*K_au), [n, N]); X = X ./ mean(X); X = X * L_au;
Y = ncx2rnd(2, sqrt(2*K_ub), [n, N]); Y = Y ./ mean(Y); Y = Y * L_ub;

%% ========== Monte Carlo: EPA (lambda = 0.5) =============================
lambda = 0.5; Pa = PMat * lambda; Pb = PMat * (1 - lambda);
epsilon  = zeta * N0^2 ./ (Pa .* X + Pb .* Y);
Gamma_AB = (eff * beta * (1 - beta) * Pa .* X .* Y) ./ ...
           (eff * beta * (1 - beta + zeta) * Y * N0 + (1 - beta) * N0 + epsilon);
Pc_Sim_EPA = mean(0.5 * log2(1 + Gamma_AB) > Rt, 2);

%% ========== Monte Carlo: OPA (lambda = 1) ===============================
lambda = 1; Pa = PMat * lambda; Pb = PMat * (1 - lambda);
epsilon  = zeta * N0^2 ./ (Pa .* X + Pb .* Y);
Gamma_AB = (eff * beta * (1 - beta) * Pa .* X .* Y) ./ ...
           (eff * beta * (1 - beta + zeta) * Y * N0 + (1 - beta) * N0 + epsilon);
Pc_Sim_OPA = mean(0.5 * log2(1 + Gamma_AB) > Rt, 2);

%% ========== Closed-Form Approximation (EPA) =============================
D_ord = 20; R_ord = 20;
Pc_Approx_EPA = zeros(1, n);
for i = 1:n
    lambda = 0.5; Pa = PVec(i) * lambda;
    A = (1 - beta + zeta) * N0 * delta_t / (1 - beta) / Pa / L_au;
    B = N0 * delta_t / eff / L_au / beta / L_ub;
    a1 = 2*K_au; a2 = 2*(1+K_au)*A; a3 = 2*(1+K_au)*B;
    a4 = K_ub + 1; a5 = sqrt(K_ub * (1+K_ub));
    Func = 0;
    for d = 0:D_ord
        for u = 0:d
            for v = 0:u
                for r = 0:R_ord
                    Func = Func + (gamma(d+D_ord)*D_ord^(1-2*d)*a1^d*a2^v*a3^(u-v) * ...
                        R_ord^(1-2*r)*nchoosek(u,v)*gamma(R_ord+r)*a5^(2*r)*2 * ...
                        (a3/2/a4)^((r-u+v+1)/2) * besselk(r-u+v+1, 2*sqrt(a3*a4/2))) / ...
                        (gamma(r+1)^2 * gamma(D_ord-d+1) * factorial(d) * factorial(u) * ...
                         gamma(R_ord-r+1) * 2^(d+u) * exp(0.5*(a1+a2)));
                end
            end
        end
    end
    Pc_Approx_EPA(i) = (1 + K_ub) * exp(-K_ub) * Func;
end

%% ========== Closed-Form Approximation (OPA) =============================
Pc_Approx_OPA = zeros(1, n);
for i = 1:n
    lambda = 1; Pa = PVec(i) * lambda;
    A = (1 - beta + zeta) * N0 * delta_t / (1 - beta) / Pa / L_au;
    B = N0 * delta_t / eff / L_au / beta / L_ub;
    a1 = 2*K_au; a2 = 2*(1+K_au)*A; a3 = 2*(1+K_au)*B;
    a4 = K_ub + 1; a5 = sqrt(K_ub * (1+K_ub));
    Func = 0;
    for d = 0:D_ord
        for u = 0:d
            for v = 0:u
                for r = 0:R_ord
                    Func = Func + ((1+K_ub)*exp(-K_ub)*gamma(d+D_ord)*D_ord^(1-2*d)* ...
                        a1^d*a2^v*a3^(u-v)*R_ord^(1-2*r)*nchoosek(u,v)*gamma(R_ord+r)* ...
                        a5^(2*r)*2*(a3/2/a4)^((r-u+v+1)/2) * ...
                        besselk(r-u+v+1,2*sqrt(a3*a4/2))) / ...
                        (gamma(r+1)^2*gamma(D_ord-d+1)*factorial(d)*factorial(u)* ...
                         gamma(R_ord-r+1)*2^(d+u)*exp(0.5*(a1+a2)));
                end
            end
        end
    end
    Pc_Approx_OPA(i) = (1 + K_ub) * exp(-K_ub) * Func;
end

%% ========== Plot ========================================================
figure('Name', 'Connection Probability: OPA vs EPA');
plot(P_dBW, Pc_Sim_OPA,   'sb',  'MarkerSize', 6, 'DisplayName', 'Simulation - OPA');
hold on;
plot(P_dBW, Pc_Approx_OPA, '-b', 'LineWidth', 1.5, 'DisplayName', 'Approximate - OPA');
plot(P_dBW, Pc_Sim_EPA,   'or',  'MarkerSize', 6, 'DisplayName', 'Simulation - EPA');
plot(P_dBW, Pc_Approx_EPA, '--r','LineWidth', 1.5, 'DisplayName', 'Approximate - EPA');

xlabel('Network Transmit Power P_a [dBW]', 'FontSize', 14, 'FontName', 'Times New Roman');
ylabel('Connection Probability',            'FontSize', 14, 'FontName', 'Times New Roman');
legend('show', 'Location', 'best', 'FontSize', 11);
grid on; box on;

%% ========== Save ========================================================
results_dir = fullfile(fileparts(mfilename('fullpath')), '..', '..', 'results', 'figures');
if ~exist(results_dir, 'dir'); mkdir(results_dir); end
savefig(fullfile(results_dir, 'CP_OPA_vs_EPA.fig'));
print( fullfile(results_dir, 'CP_OPA_vs_EPA.png'), '-dpng', '-r300');

function [L, K, alpha_eff] = compute_air_to_ground_channel(d_horiz, H, alpha_L, alpha_N, a, b, km, kM)
% COMPUTE_AIR_TO_GROUND_CHANNEL  Computes large-scale path loss and Rician K-factor
%   for an Air-to-Ground (A2G) link based on elevation angle statistics.
%
%   [L, K, alpha_eff] = compute_air_to_ground_channel(d_horiz, H, alpha_L, alpha_N, a, b, km, kM)
%
%   The function uses the probabilistic LoS model for urban environments to
%   compute the expected path loss and the elevation-angle-dependent Rician
%   K-factor (in linear scale) for the A2G channel.
%
%   Inputs:
%       d_horiz  - Horizontal distance between ground node and UAV (m)
%       H        - UAV altitude (m)
%       alpha_L  - Path-loss exponent for LoS links (typical: 2.0)
%       alpha_N  - Path-loss exponent for NLoS links (typical: 3.5-4.0)
%       a        - LoS probability parameter (urban: 0.28)
%       b        - LoS probability parameter (urban: 9.61)
%       km       - Minimum Rician K-factor in dB (e.g., 1 dB)
%       kM       - Maximum Rician K-factor in dB (e.g., 10 dB)
%
%   Outputs:
%       L        - Average large-scale path loss (linear, dimensionless)
%       K        - Rician K-factor (linear scale, elevation-dependent)
%       alpha_eff- Effective path-loss exponent (weighted LoS/NLoS)
%
%   Model:
%       - Elevation angle: theta = asin(H / d_3D)
%       - LoS probability: P_LoS = 1 / (1 + a * exp(-b * (theta - a)))
%       - Effective exponent: alpha_eff = (alpha_L - alpha_N) * P_LoS + alpha_N
%       - Path loss: L = d_3D^(-alpha_eff)
%       - K-factor [dB]: K_dB = km + (kM - km) * (2/pi) * theta
%
%   Reference: See system model in the published paper.
%
%   See also: generate_rician_channel, generate_rayleigh_channel

    d_3D    = sqrt(d_horiz^2 + H^2);           % 3D Euclidean distance
    theta   = asin(H / d_3D);                  % Elevation angle (radians)

    % LoS probability (Sigmoid-based urban model)
    P_los   = 1 / (1 + a * exp(-b * (theta - a)));

    % Effective path-loss exponent
    alpha_eff = (alpha_L - alpha_N) * P_los + alpha_N;

    % Average large-scale path loss
    L = d_3D^(-alpha_eff);

    % Rician K-factor (linear) — elevation-angle dependent
    K_dB = km + (kM - km) * (2 / pi) * theta;
    K    = db2pow(K_dB);
end


function [X, lambda_X] = generate_rician_channel(K, L, n_rows, n_cols)
% GENERATE_RICIAN_CHANNEL  Generates normalized Rician fading channel realizations.
%
%   [X, lambda_X] = generate_rician_channel(K, L, n_rows, n_cols)
%
%   Generates a matrix of channel power gains |h|^2 for a Rician-fading
%   Air-to-Ground link, normalized so that E[X] = L (the large-scale path loss).
%
%   The small-scale fading follows a non-central chi-squared distribution
%   (2 degrees of freedom) with non-centrality parameter lambda_X = sqrt(2*K).
%
%   Inputs:
%       K       - Rician K-factor (linear)
%       L       - Large-scale path loss (average power, linear)
%       n_rows  - Number of rows in output (e.g., simulation steps)
%       n_cols  - Number of Monte Carlo samples per step
%
%   Outputs:
%       X       - Matrix of channel power gains (n_rows x n_cols), scaled to E[X] = L
%       lambda_X- Non-centrality parameter used in the ncx2 distribution
%
%   Model:
%       lambda_X = sqrt(2*K)
%       X_raw ~ ncx2(2, lambda_X)
%       X = (X_raw / mean(X_raw)) * L    [normalized so E[X] = L]
%
%   See also: generate_rayleigh_channel, compute_air_to_ground_channel

    lambda_X = sqrt(2 * K);
    X_raw    = ncx2rnd(2, lambda_X, [n_rows, n_cols]);
    X        = (X_raw ./ mean(X_raw(:))) * L;
end


function V = generate_rayleigh_channel(L, n_rows, n_cols)
% GENERATE_RAYLEIGH_CHANNEL  Generates exponential (Rayleigh-power) channel realizations.
%
%   V = generate_rayleigh_channel(L, n_rows, n_cols)
%
%   Generates a matrix of channel power gains |h|^2 for a Rayleigh-fading
%   Ground-to-Ground (G2G) link.  The power gain follows an exponential
%   distribution with mean L (the large-scale path loss).
%
%   Inputs:
%       L       - Large-scale path loss (average power, linear)
%       n_rows  - Number of rows in output
%       n_cols  - Number of Monte Carlo samples per row
%
%   Output:
%       V       - Matrix of channel power gains (n_rows x n_cols), E[V] = L
%
%   See also: generate_rician_channel

    V = L * exprnd(1, [n_rows, n_cols]);
end

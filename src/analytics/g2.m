function y = g2(lambda, b, R)
% G2  Approximation of E[log(X + b)] for a non-central chi-squared RV X ~ ncx2(2, lambda).
%
%   y = g2(lambda, b, R) computes the R-term series approximation of
%   E[log(X + b)] where X is a normalized noncentral chi-squared RV with
%   2 degrees of freedom and non-centrality parameter lambda.
%
%   This function is used in the ESR lower-bound approximation to account
%   for the additive noise offset b in the denominator of the SNR expression.
%
%   Inputs:
%       lambda  - Non-centrality parameter (lambda >= 0)
%       b       - Positive scalar offset: (1 - beta) / (eta * beta * (1 - beta + zeta))
%                 (i.e., the ratio capturing noise and SWIPT parameters)
%       R       - Number of series terms (R = 5..30 is typically sufficient)
%
%   Output:
%       y       - Approximation of E[log(X + b)]
%
%   Note: Calls PHI(r, b) internally.
%
%   Reference: See Eq. (XX) in the published paper.
%
%   See also: g1, PHI, PSI

    Temp = 0;
    for r = 0:R
        Temp = Temp + (gamma(R + r) * R^(1 - 2*r) * (lambda/4)^r ...
                       / (gamma(r + 1)^2) / gamma(R - r + 1)) * PHI(r, b);
    end
    y = 0.5 * exp(-lambda/2) * Temp;
end

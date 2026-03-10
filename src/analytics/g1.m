function y = g1(lambda, R)
% G1  Approximation of E[log(X)] for a non-central chi-squared RV X ~ ncx2(2, lambda).
%
%   y = g1(lambda, R) computes the R-term series approximation of E[log(X)]
%   where X follows a noncentral chi-squared distribution with 2 degrees of
%   freedom and non-centrality parameter lambda.
%
%   This result is used in the lower-bound approximation of the ergodic
%   secrecy rate (ESR) for the A->U and U->B links with Rician fading.
%
%   Inputs:
%       lambda  - Non-centrality parameter of the ncx2 distribution (lambda >= 0)
%       R       - Number of terms in the series expansion (controls accuracy;
%                 R = 5..30 is typically sufficient)
%
%   Output:
%       y       - Approximation of E[log(X)]
%
%   Reference: See Eq. (XX) in the published paper.
%
%   See also: g2, PHI, PSI

    Temp = 0;
    for r = 0:R
        Temp = Temp + gamma(R + r) * R^(1 - 2*r) * (lambda/2)^r ...
                      * (psi(r + 1) + log(2)) ...
                      / gamma(r + 1) / gamma(R - r + 1);
    end
    y = exp(-lambda/2) * Temp;
end

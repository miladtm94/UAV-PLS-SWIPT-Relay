function y = g3(lambda, R)
% G3  Series approximation of E[1/X] for a non-central chi-squared RV X ~ ncx2(2, lambda).
%
%   y = g3(lambda, R) computes the R-term series approximation of E[1/X]
%   where X is a normalized noncentral chi-squared RV with 2 degrees of
%   freedom and non-centrality parameter lambda.
%
%   Inputs:
%       lambda  - Non-centrality parameter (lambda >= 0)
%       R       - Number of series terms (default 50 if not supplied)
%
%   Output:
%       y       - Approximation of E[1/X]
%
%   Reference: See Eq. (XX) in the published paper.
%
%   See also: g1, g2

    if nargin < 2
        R = 50;
    end

    y = 0;
    for r = 0:R
        y = y + gamma(R + r) * gamma(r) * R^(1 - 2*r) * (lambda/2)^r ...
              * exp(-lambda/2) / gamma(R - r + 1) / gamma(r + 1)^2 / 2;
    end
end

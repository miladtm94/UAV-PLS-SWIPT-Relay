function y = PSI(m, r, beta)
% PSI  Computes the PSI(m, r, beta) auxiliary series used in the closed-form ESR.
%
%   y = PSI(m, r, beta) evaluates the PSI summation that appears in the
%   second-order moment integral of the U->B channel SNR expression.
%
%   Inputs:
%       m     - Lower summation index (non-negative integer)
%       r     - Upper summation index (non-negative integer, r >= m)
%       beta  - Power splitting ratio (0 < beta < 1)
%
%   Output:
%       y     - Scalar value of PSI(m, r, beta)
%
%   Reference: See Eq. (XX) in the published paper.
%
%   See also: g2, PHI

    y = 0;
    for k = 1:(r - m)
        y = y + factorial(k - 1) * (-beta/2)^(r - m - k);
    end
end

function y = PHI(r, b)
% PHI  Computes the PHI(r, b) auxiliary function used in the ESR closed-form expression.
%
%   y = PHI(r, b) evaluates a summation involving the incomplete gamma
%   function and the Meijer G-function, arising in the closed-form
%   derivation of E[log(X + b)] for a non-central chi-squared RV.
%
%   Inputs:
%       r   - Non-negative integer order
%       b   - Positive scalar: the additive offset (ratio of SWIPT parameters)
%
%   Output:
%       y   - Scalar value of PHI(r, b)
%
%   Note: Requires MATLAB Symbolic Math Toolbox for meijerG().
%
%   Reference: See Eq. (XX) in the published paper.
%
%   See also: g2, PSI

    temp = 0;
    for s = 0:r
        temp = temp + nchoosek(r, s) * (-b)^(r - s) ...
                    * (2^(1 + s) * (log(b) * gamma(s + 1) ...
                       * gammainc(b/2, s + 1, 'upper') ...
                       + meijerG([], [1, 1], [0, 0, 1 + s], [], b/2)));
    end
    y = exp(b/2) * temp;
end

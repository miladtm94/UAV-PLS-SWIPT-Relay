function y = I1(t, beta)
% I1  Evaluates the I1(t, beta) auxiliary integral used in the ESR derivation.
%
%   y = I1(t, beta) computes a closed-form expression for the integral
%   that arises in the ergodic capacity analysis of the U->B Rician channel
%   with a power-splitting SWIPT receiver.
%
%   Inputs:
%       t     - Non-negative integer (order parameter)
%       beta  - Power splitting ratio at the UAV (0 < beta < 1)
%
%   Output:
%       y     - Scalar value of I1(t, beta)
%
%   Reference: See Eq. (XX) in the published paper.
%
%   See also: g2, Ei, PHI

    y = 0;
    if t == 0
        y = -2 * exp(-beta/2) * Ei(-beta/2);
    end
    if t > 0
        Temp = 0;
        for m = 0:(t - 1)
            Temp = Temp + (factorial(t) / factorial(t - m)) ...
                        * (((-1)^(t - m - 1) / ((2/beta)^(t - m))) ...
                           * exp(-beta/2) * Ei(-beta/2));
            for k = 1:(t - m)
                Temp = Temp + factorial(k - 1) / ((-2/beta)^(t - m - k));
            end
        end
        y = 2^(t + 1) * (Temp + factorial(t) * exp(beta/2) * Ei(-beta/2));
    end
end

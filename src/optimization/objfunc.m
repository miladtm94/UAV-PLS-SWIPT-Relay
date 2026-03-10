function f = objfunc(beta, coef_a, coef_b)
% OBJFUNC  Polynomial ratio objective for PSR (beta) optimization.
%
%   f = objfunc(beta, coef_a, coef_b) evaluates the secrecy-rate expression
%   written as a ratio of two polynomials in beta, enabling optimization of
%   the Power Splitting Ratio (PSR) beta in [0, 1].
%
%   The secrecy rate is:
%       Rs(beta) = 0.5 * log2[ Num(beta) / Den(beta) ]
%   where:
%       Num(beta) = a(1) + a(2)*beta + a(3)*beta^2 + a(4)*beta^3 + ...
%       Den(beta) = b(1) + b(2)*beta + b(3)*beta^2 + b(4)*beta^3 + b(5)*beta^4 + ...
%
%   The polynomial coefficients coef_a and coef_b are computed symbolically
%   from the SNR expressions (see Ultimate_Simulation.m for derivation).
%
%   Inputs:
%       beta   - Power splitting ratio scalar (or array), 0 < beta < 1
%       coef_a - Numerator polynomial coefficients [a1, a2, a3, a4, ...]
%       coef_b - Denominator polynomial coefficients [b1, b2, b3, b4, b5, ...]
%
%   Output:
%       f      - Ratio Num(beta)/Den(beta) (NOT yet in log scale; maximizing f
%                is equivalent to maximizing Rs)
%
%   Note: This function is used with MATLAB's optimproblem / solve framework.
%
%   Reference: See Eq. (XX) and the Dinkelbach-method discussion in the paper.
%
%   See also: obj, Ultimate_Simulation

    Num = @(b) coef_a(1) + coef_a(2)*b + coef_a(3)*b.^2 + coef_a(4)*b.^3;
    Den = @(b) coef_b(1) + coef_b(2)*b + coef_b(3)*b.^2 ...
                         + coef_b(4)*b.^3 + coef_b(5)*b.^4;
    f = Num(beta) ./ Den(beta);
end

function out = coeff(R, r)
% COEFF  Computes the series expansion coefficient used in connection probability.
%
%   out = coeff(R, r) returns the coefficient:
%       gamma(R + r) * R^(1 - 2*r) / (gamma(r+1))^2 / gamma(R - r + 1)
%
%   This coefficient appears in the closed-form Marcum-Q function series
%   expansion used to derive the analytical connection probability.
%
%   Inputs:
%       R   - Series truncation order (positive integer)
%       r   - Current summation index (0 <= r <= R)
%
%   Output:
%       out - Scalar coefficient value
%
%   Reference: See Eq. (XX) in the published paper.
%
%   See also: g1, g2

    out = gamma(R + r) * R^(1 - 2*r) / (gamma(r + 1))^2 / gamma(R - r + 1);
end

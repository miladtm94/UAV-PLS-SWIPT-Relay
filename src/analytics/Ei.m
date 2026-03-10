function y = Ei(x)
% EI  Computes the exponential integral Ei(x) = -E1(-x).
%
%   y = Ei(x) returns the exponential integral for real x > 0, defined as:
%       Ei(x) = -PV integral_{-x}^{inf} exp(-t)/t dt
%            = -expint(-x)
%
%   This function wraps MATLAB's built-in expint() to provide the standard
%   Ei(x) sign convention used in the analytical expressions for the
%   ergodic secrecy rate.
%
%   Input:
%       x   - Real scalar or array (x > 0 for the standard definition)
%
%   Output:
%       y   - Exponential integral Ei(x)
%
%   Reference: Used in the closed-form ESR lower bound; see the paper.
%
%   See also: expint, g2

    y = -expint(-x);
end

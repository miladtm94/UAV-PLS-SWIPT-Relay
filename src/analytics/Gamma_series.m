function temp = Gamma_series(s, b)
% GAMMA_SERIES  Evaluates a truncated series representation of a Gamma-related function.
%
%   temp = Gamma_series(s, b) computes:
%       sum_{m=0}^{s} s! * exp(-b/2) * (b/2)^m / m!
%
%   This helper function appears in the evaluation of certain integral
%   expressions in the connection probability analysis.
%
%   Inputs:
%       s   - Non-negative integer (upper summation limit)
%       b   - Positive real scalar
%
%   Output:
%       temp - Scalar result of the series
%
%   Note: This function was originally named Gamma.m but has been renamed
%   to Gamma_series.m to avoid shadowing MATLAB's built-in gamma() function.
%
%   See also: gamma, coeff

    temp = 0;
    for m = 0:s
        temp = temp + factorial(s) * exp(-b/2) * ((b/2)^m) / factorial(m);
    end
end

function y = obj(c1, c2, c3, c4, c5, lambda)
% OBJ  Instantaneous secrecy rate objective function for power allocation optimization.
%
%   y = obj(c1, c2, c3, c4, c5, lambda) evaluates the (negated) instantaneous
%   secrecy rate as a function of the power allocation factor lambda, for use
%   with MATLAB's fminunc/fmincon optimizers.
%
%   The secrecy rate is expressed as:
%       Rs(lambda) = 0.5 * log2[ (1 + c1*lambda) / (1 + c2*lambda/(1+c3*(1-lambda))
%                                                       + c4*lambda/(1+c5*(1-lambda))) ]
%
%   where the composite constants c1..c5 capture channel gains, power, and
%   SWIPT parameters (see the paper for full definitions).
%
%   The function returns the NEGATIVE of Rs so that minimizing y maximizes Rs.
%
%   Inputs:
%       c1   - Effective desired-link SNR coefficient (scalar or array)
%       c2   - Eavesdropper direct-path SNR coefficient (lambda-related)
%       c3   - Eavesdropper direct-path jamming coefficient
%       c4   - Eavesdropper relay-path SNR coefficient
%       c5   - Eavesdropper relay-path jamming coefficient
%       lambda - Power allocation factor [0, 1] (fraction allocated to Alice->U link)
%
%   Output:
%       y    - Negative instantaneous secrecy rate (for minimization)
%
%   Reference: See Eq. (XX) in the published paper.
%
%   See also: optimal_power_allocation, objfunc

    y = -((1 + c1 .* lambda) ./ ...
          (1 + c2 .* lambda ./ (1 + c3 .* (1 - lambda)) ...
             + c4 .* lambda ./ (1 + c5 .* (1 - lambda))));
end

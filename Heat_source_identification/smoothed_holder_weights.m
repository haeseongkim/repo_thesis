% Function to calculate smoothed Hölder weights
% This function computes weights based on the smoothed Hölder approach,
% which incorporates a smoothing term (epsilon) and a power parameter (p).

function z = smoothed_holder_weights(x, epsilon, p)
    % Input parameters:
    % x - The input array or matrix for which weights are to be calculated.
    % epsilon - A small positive value to ensure stability and smoothness.
    % p - The power parameter, influencing the weight calculation.

    % Compute the weights.
    % The operation is element-wise, denoted by .^, which ensures that the
    % computation is performed on each element of the array x individually.
    % The formula used is: (x^2 + epsilon^2)^(p/2 - 1),
    % where x^2 and epsilon^2 ensure non-negativity,
    % and (p/2 - 1) adjusts the rate of change of weights with respect to x.
    z = (x.^2 + epsilon^2).^(p/2 - 1);

    % The output z is an array of the same size as x,
    % containing the calculated weights for each element of x.
end

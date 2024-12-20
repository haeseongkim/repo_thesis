% discontinuity adaptive Markov Random Field: line process
% input eta: absolute difference between adjacent pixel, threshold: magnitude
function h = MRF_line_process(eta, threshold)
    % h_line_process returns 1 where eta below threshold, and 0 otherwise
    % Vectorized approach to handle arrays or single values of eta without using if statements

    % Perform the logical comparison and convert the result to double
    % MATLAB treats true as 1 and false as 0 automatically
    h = double(abs(eta) < threshold); % 
    
end

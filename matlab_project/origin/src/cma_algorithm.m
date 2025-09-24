function [y, w, error_curve] = cma_algorithm(x, step_size, num_taps, num_iterations)
% CMA_ALGORITHM - Constant Modulus Algorithm Implementation
%
% The Constant Modulus Algorithm (CMA) is a blind equalization technique
% used to remove intersymbol interference from communication signals.
%
% Inputs:
%   x           - Input signal (received signal)
%   step_size   - Algorithm step size (mu)
%   num_taps    - Number of equalizer taps
%   num_iterations - Number of iterations to run
%
% Outputs:
%   y           - Equalized output signal
%   w           - Final equalizer weights
%   error_curve - CMA error over iterations
%
% Author: DSP Course Design
% Date: 2025

    % Initialize parameters
    N = length(x);
    
    % Initialize equalizer weights (centered tap initialization)
    w = zeros(num_taps, 1);
    center_tap = ceil(num_taps/2);
    w(center_tap) = 1;  % Center tap initialization
    
    % Output and error storage
    y = zeros(N, 1);
    error_curve = zeros(num_iterations, 1);
    
    % CMA constant (for QPSK, R2 = 1)
    R2 = 1;
    
    % Main CMA loop
    for iter = 1:num_iterations
        total_error = 0;
        
        for n = num_taps:N
            % Extract input vector
            x_vec = x(n:-1:n-num_taps+1);
            
            % Equalizer output
            y(n) = w' * x_vec;
            
            % CMA error
            error = y(n) * (R2 - abs(y(n))^2);
            total_error = total_error + abs(error)^2;
            
            % Weight update
            w = w + step_size * conj(error) * x_vec;
        end
        
        % Store average error for this iteration
        error_curve(iter) = total_error / (N - num_taps + 1);
        
        % Optional: Early stopping if converged
        if iter > 1 && abs(error_curve(iter) - error_curve(iter-1)) < 1e-8
            error_curve = error_curve(1:iter);
            break;
        end
    end
end
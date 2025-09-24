function [X] = fft_implementation(x)
% FFT_IMPLEMENTATION - Radix-2 Decimation-In-Time FFT Algorithm
% 
% This function implements the Fast Fourier Transform using the
% radix-2 decimation-in-time algorithm.
%
% Input:
%   x - Input signal (must be power of 2 length)
%
% Output:
%   X - FFT of input signal
%
% Author: DSP Course Design
% Date: 2025

    N = length(x);
    
    % Check if N is power of 2
    if mod(log2(N), 1) ~= 0
        error('Input length must be a power of 2');
    end
    
    % Bit-reversal permutation
    x = bit_reverse(x, N);
    
    % Number of stages
    stages = log2(N);
    
    % FFT computation
    for stage = 1:stages
        M = 2^stage;          % Size of each DFT
        M_half = M/2;
        
        % Twiddle factor
        W = exp(-1j * 2 * pi / M);
        
        for k = 0:(N/M-1)
            for j = 0:(M_half-1)
                idx1 = k*M + j + 1;
                idx2 = idx1 + M_half;
                
                % Butterfly computation
                temp = x(idx2) * W^j;
                x(idx2) = x(idx1) - temp;
                x(idx1) = x(idx1) + temp;
            end
        end
    end
    
    X = x;
end

function x_reversed = bit_reverse(x, N)
% BIT_REVERSE - Perform bit-reversal permutation
    
    x_reversed = zeros(size(x));
    bits = log2(N);
    
    for i = 0:N-1
        reversed_i = 0;
        temp_i = i;
        
        for b = 0:bits-1
            reversed_i = reversed_i * 2 + mod(temp_i, 2);
            temp_i = floor(temp_i / 2);
        end
        
        x_reversed(reversed_i + 1) = x(i + 1);
    end
end
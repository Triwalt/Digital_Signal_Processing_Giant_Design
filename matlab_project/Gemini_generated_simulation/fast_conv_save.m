function y = fast_conv_save(x, h, N)
% fast_conv_save: Implements fast convolution using the Overlap-Save method.
% This method is also efficient for convolving a long signal 'x' with 'h'.
%
% Arguments:
%   x : The long input signal vector.
%   h : The filter impulse response vector.
%   N : The FFT size (block size), must be a power of 2 and > length(h)-1.

    M = length(h);
    Nx = length(x);
    
    % N must be large enough to contain the filter response.
    if N <= M-1
        error('FFT size N must be greater than filter length M-1.');
    end
    
    % The number of valid output points from each processed block.
    L = N - M + 1;
    
    % Pad the filter 'h' to the FFT size and compute its FFT once.
    h_padded = [h, zeros(1, N - M)];
    H = my_fft(h_padded);

    % --- BUG FIX ---
    % The previous logic for padding and block calculation was incorrect,
    % leading to an output vector that was too short. The corrected logic
    % robustly calculates the required padding and number of blocks.
    
    % 1. Pre-pend M-1 zeros to the input signal to handle initial overlap.
    x_prepended = [zeros(1, M - 1), x];
    
    % 2. Calculate the number of blocks needed to cover the entire prepended signal.
    %    Each block produces L valid output samples.
    num_blocks = ceil(length(x_prepended) / L);
    
    % 3. Pad the signal at the end so its total length is suitable for block processing.
    target_len = (num_blocks - 1) * L + N;
    x_padded = [x_prepended, zeros(1, target_len - length(x_prepended))];
    
    % Initialize the output vector to hold all valid samples.
    y_full = zeros(1, num_blocks * L);
    
    % Process each block.
    for i = 1:num_blocks
        % Extract the current overlapping block of size N.
        start_idx = (i-1)*L + 1;
        end_idx = start_idx + N - 1;
        block_x = x_padded(start_idx:end_idx);
        
        % Compute the FFT of the block.
        Block_X = my_fft(block_x);
        
        % Perform frequency-domain multiplication.
        Block_Y = Block_X .* H;
        
        % Compute the inverse FFT.
        block_y = my_ifft(Block_Y);
        
        % "Save" step: Discard the first M-1 points and keep the valid L points.
        valid_y = block_y(M:N);
        
        % Place the valid samples into the correct position in the output vector.
        y_start = (i-1)*L + 1;
        y_end = i*L;
        y_full(y_start:y_end) = valid_y;
    end
    
    % Trim the fully constructed output vector to the correct final length (Nx + M - 1).
    % The y_full vector is now guaranteed to be long enough.
    y = y_full(1:(Nx + M - 1));
end


function y = fast_conv_add(x, h, L)
% fast_conv_add: Implements fast convolution using the Overlap-Add method.
% This method is efficient for convolving a very long signal 'x' with a 
% shorter filter 'h'.
%
% Arguments:
%   x : The long input signal vector.
%   h : The filter impulse response vector.
%   L : The processing block length for the input signal 'x'.

    M = length(h);
    Nx = length(x);
    
    % Determine the FFT size. It must be large enough to avoid time-domain aliasing.
    % For my_fft, we pad this to the next power of 2.
    N_fft_min = L + M - 1;
    N = 2^nextpow2(N_fft_min);
    
    % Pad the filter 'h' to the FFT size and compute its FFT once.
    h_padded = [h, zeros(1, N - M)];
    H = my_fft(h_padded);
    
    % Determine the number of blocks needed to process the entire signal.
    num_blocks = ceil(Nx / L);
    
    % Pad the input signal 'x' so its length is a multiple of L.
    x_padded = [x, zeros(1, L * num_blocks - Nx)];
    
    % Initialize the output vector and an overlap buffer.
    y = zeros(1, Nx + M - 1);
    
    % Process each block of the input signal.
    for i = 1:num_blocks
        % Extract the current block from the padded input signal.
        start_idx = (i-1)*L + 1;
        end_idx = i*L;
        block_x = x_padded(start_idx:end_idx);
        
        % Pad the block to the FFT size and compute its FFT.
        block_x_padded = [block_x, zeros(1, N - L)];
        Block_X = my_fft(block_x_padded);
        
        % Perform frequency-domain multiplication.
        Block_Y = Block_X .* H;
        
        % Compute the inverse FFT to get the block convolution result.
        block_y = my_ifft(Block_Y);
        
        % Perform the Overlap-Add operation.
        y_start = (i-1)*L + 1;
        y_end = y_start + N - 1;
        
        % Add the convoluted block to the corresponding section of the output vector.
        % This correctly handles the overlapping regions between blocks.
        if y_end > length(y)
            y_end_limited = length(y);
            y(y_start:y_end_limited) = y(y_start:y_end_limited) + block_y(1:(y_end_limited - y_start + 1));
        else
            y(y_start:y_end) = y(y_start:y_end) + block_y;
        end
    end
    
    % Trim the output vector to the correct final length (Nx + M - 1).
    y = y(1:(Nx + M - 1));
end

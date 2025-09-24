function Y = my_fft(x)
% my_fft: Implements a recursive Radix-2 DIT FFT algorithm.
% It pads the input to the next power of 2 to handle arbitrary length inputs,
% as required by the basic radix-2 algorithm.

    N = length(x);
    
    % Calculate the next power of 2 for padding
    N_pow2 = 2^nextpow2(N);
    if N ~= N_pow2
        % Pad with zeros if the length is not a power of 2
        x = [x, zeros(1, N_pow2 - N)];
        N = N_pow2;
    end

    % Base case of the recursion: a 1-point FFT is the sample itself
    if N == 1
        Y = x;
        return;
    end

    % Recursive step: split into even and odd parts
    x_even = x(1:2:N);
    x_odd = x(2:2:N);

    % Perform FFT on each half
    G = my_fft(x_even);
    H = my_fft(x_odd);

    % Manually calculate twiddle factors, W_N^k
    % W_N^k = exp(-1j * 2 * pi * k / N)
    % This loop implements the calculation based on the formula.
    k = 0:(N/2 - 1);
    W = exp(-1j * 2 * pi * k / N);
    
    % Combine the results using the butterfly operation
    WH = W .* H;
    Y = [G + WH, G - WH];
end

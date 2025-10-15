function y = my_ifft(X)
% my_ifft: Implements the Inverse Fast Fourier Transform.
% It uses the custom my_fft function based on the mathematical property:
% IFFT(X) = (1/N) * conj(FFT(conj(X)))
% This is an efficient way to implement IFFT when an FFT function is available.

    N = length(X);

    % Enforce power-of-two length to avoid unintended padding inside my_fft
    if N ~= 2^nextpow2(N)
        error('my_ifft: input length must be a power of 2. Got N = %d.', N);
    end
    
    % Step 1: Compute the complex conjugate of the input sequence.
    conj_X = conj(X);
    
    % Step 2: Apply the forward FFT algorithm to the conjugated sequence.
    fft_of_conj = my_fft(conj_X);
    
    % Step 3: Compute the complex conjugate of the result and scale by 1/N.
    y = conj(fft_of_conj) / N;
end

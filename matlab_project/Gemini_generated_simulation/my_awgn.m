function y = my_awgn(x, snr_db, power_mode)
% my_awgn: Adds white Gaussian noise to signal to achieve specified SNR
% 
% INPUTS:
%   x         - Input signal
%   snr_db    - Signal-to-noise ratio in dB
%   power_mode - 'measured' to measure signal power, or numeric value
%
% OUTPUT:
%   y         - Signal with added noise
%
% ALGORITHM:
% 1. Calculate signal power based on power_mode
% 2. Determine noise power from SNR requirement
% 3. Generate complex white Gaussian noise
% 4. Scale noise to achieve target SNR
% 5. Add noise to signal

    % Ensure input is a row vector for consistency
    if size(x, 1) > 1
        x = x.';
        transpose_output = true;
    else
        transpose_output = false;
    end
    
    % Calculate signal power
    if strcmp(power_mode, 'measured')
        % Measure the actual signal power
        signal_power = mean(abs(x).^2);
    else
        % Use specified signal power (assume normalized to 1 if not specified)
        signal_power = 1;
    end
    
    % Convert SNR from dB to linear scale
    snr_linear = 10^(snr_db / 10);
    
    % Calculate required noise power
    noise_power = signal_power / snr_linear;
    
    % Generate complex white Gaussian noise
    % For complex signals: noise has power split equally between real and imaginary parts
    if ~isreal(x)
        % Complex noise: variance = noise_power/2 for each component
        noise_real = sqrt(noise_power/2) * randn(size(x));
        noise_imag = sqrt(noise_power/2) * randn(size(x));
        noise = noise_real + 1j * noise_imag;
    else
        % Real noise: variance = noise_power
        noise = sqrt(noise_power) * randn(size(x));
    end
    
    % Add noise to signal
    y = x + noise;
    
    % Restore original orientation if needed
    if transpose_output
        y = y.';
    end
    
    % Verification (optional debug output)
    % actual_snr = 10*log10(mean(abs(x).^2) / mean(abs(noise).^2));
    % fprintf('Target SNR: %.2f dB, Actual SNR: %.2f dB\n', snr_db, actual_snr);
end

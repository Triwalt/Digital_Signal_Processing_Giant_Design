%% Main Simulation Script for Optical Communication Signal Processing
%
% This script demonstrates and verifies the custom-built modules:
% 1. Verifies my_fft against MATLAB's built-in fft.
% 2. Verifies fast_conv_add and fast_conv_save against MATLAB's conv.
% 3. Simulates a 4-QAM communication link with an ISI channel and uses a
%    Constant Modulus Algorithm (CMA) equalizer to mitigate the distortion.

clear; clc; close all;

%% Part 1: Verification of my_fft and my_ifft

fprintf('--- Part 1: Verifying my_fft and my_ifft ---\n');

% Generate a random complex signal for testing.
N_test = 100; % Use an arbitrary length to test padding.
test_signal = randn(1, N_test) + 1j * randn(1, N_test);

% Run custom FFT. Our function will pad it to N_fft = 128.
N_fft = 2^nextpow2(N_test);
fft_custom = my_fft(test_signal);

% Run MATLAB's built-in FFT for comparison.
fft_matlab = fft(test_signal, N_fft);

% Calculate and display the maximum absolute error.
fft_error = max(abs(fft_custom - fft_matlab));
fprintf('Max error between my_fft and MATLAB fft: %e\n', fft_error);

% Verify IFFT by performing a round trip (FFT -> IFFT).
ifft_custom = my_ifft(fft_custom);
ifft_error = max(abs(ifft_custom(1:N_test) - test_signal));
fprintf('Max error for my_ifft (round trip): %e\n\n', ifft_error);


%% Part 2: Verification of Fast Convolution Methods

fprintf('--- Part 2: Verifying Fast Convolution ---\n');

% Generate a long signal and a shorter filter for testing.
x_long = randn(1, 1000); 
h_short = randn(1, 50);

% 1. Use standard time-domain convolution as the ground truth.
y_ref = conv(x_long, h_short);

% 2. Test the Overlap-Add method.
L = 128; % A suitable block size for processing.
y_add = fast_conv_add(x_long, h_short, L);
conv_add_error = max(abs(y_add - y_ref));
fprintf('Max error for Overlap-Add vs conv: %e\n', conv_add_error);

% 3. Test the Overlap-Save method.
N_save = 256; % A suitable FFT size for overlap-save.
y_save = fast_conv_save(x_long, h_short, N_save);
conv_save_error = max(abs(y_save - y_ref)); 
fprintf('Max error for Overlap-Save vs conv: %e\n\n', conv_save_error);


%% Part 3: CMA Equalizer Simulation for 4QAM

fprintf('--- Part 3: CMA Equalizer Simulation ---\n');

% --- Simulation Parameters ---
num_symbols = 20000;    % Number of symbols to transmit
M = 4;                  % Modulation order for 4-QAM
snr_db = 25;            % Signal-to-Noise Ratio in dB

% --- Channel and Equalizer Parameters ---
channel = [1, 0.4*exp(1j*pi/6), 0.1*exp(-1j*pi/4)]; % Example ISI channel
eq_taps = 21;           % Number of equalizer filter taps
mu = 0.001;             % CMA algorithm step size (learning rate)
block_size = 512;       % Processing block size

% --- 4-QAM Signal Generation ---
bits_per_symbol = log2(M);
num_bits = num_symbols * bits_per_symbol;
tx_bits = randi([0 1], 1, num_bits);
% Map bits to 4-QAM symbols { (±1 ± j) / sqrt(2) } for unit power.
s = 1/sqrt(2) * ( (1-2*tx_bits(1:2:end)) + 1j*(1-2*tx_bits(2:2:end)) );

% --- Channel Simulation ---
% Pass the signal through the ISI channel.
rx_signal_isi = conv(s, channel);
% Add Additive White Gaussian Noise (AWGN).
rx_signal = awgn(rx_signal_isi, snr_db, 'measured');

% --- CMA Equalization ---
% This implements a block-based CMA where the filter is updated once per block.
fprintf('Running Block CMA Equalizer...\n');

% Initialize equalizer weights with a center-spike.
w = zeros(1, eq_taps);
w(ceil(eq_taps/2)) = 1;
% For our unit-power 4QAM, the constant modulus radius squared R2 is 1.
R2 = 1;

num_blocks_cma = floor(length(rx_signal) / block_size);
equalized_signal = zeros(1, num_blocks_cma * block_size);

% For continuous filtering using overlap-save method
rx_padded = [zeros(1, eq_taps-1), rx_signal];
N_conv = 2^nextpow2(block_size + eq_taps - 1);

for i = 1:num_blocks_cma
    % --- FILTERING STEP using custom FFT ---
    % This block demonstrates how fast convolution can be used for filtering.
    % We get an input segment that will produce one clean output block.
    input_start = (i-1)*block_size + 1;
    input_end = input_start + block_size + eq_taps - 1 - 1;
    if input_end > length(rx_padded), continue; end
    input_segment = rx_padded(input_start:input_end);
    
    % Perform convolution in frequency domain
    w_padded = [w, zeros(1, N_conv - eq_taps)];
    W_fft = my_fft(w_padded);
    input_fft = my_fft([input_segment, zeros(1, N_conv - length(input_segment))]);
    output_fft = input_fft .* W_fft;
    y_full = my_ifft(output_fft);
    y_block = y_full(eq_taps : eq_taps + block_size - 1); % Discard aliased part

    % --- UPDATE STEP using CMA logic ---
    % The filter 'w' is updated sample-by-sample based on the output of the
    % current block. This new 'w' will be used for the next block.
    for k = 1:block_size
        % Reconstruct the filter's input vector for this sample.
        rx_idx = (i-1)*block_size + k + (eq_taps - 1);
        filter_input = rx_padded(rx_idx:-1:rx_idx - eq_taps + 1);
        
        y_k = y_block(k);
        e_k = y_k * (R2 - abs(y_k)^2); % CMA error calculation
        w = w + mu * e_k * conj(filter_input); % Update weights
    end
    
    % Store the equalized block.
    equalized_signal((i-1)*block_size+1 : i*block_size) = y_block;
end

fprintf('Equalization complete.\n');

% --- Visualization of Results ---
figure('Position', [100, 100, 1200, 400]);
sgtitle('4-QAM CMA Equalization Results');

% Plot Received Signal Constellation (distorted)
subplot(1, 3, 1);
plot(real(rx_signal), imag(rx_signal), '.b');
title('Received Signal (with ISI & Noise)');
xlabel('In-Phase'); ylabel('Quadrature');
axis([-2 2 -2 2]); grid on; axis square;

% Plot Equalized Signal Constellation
subplot(1, 3, 2);
% Discard the initial transient part of the signal before plotting.
transient_samples = 2000;
plot(real(equalized_signal(transient_samples:end)), imag(equalized_signal(transient_samples:end)), '.r');
title('Equalized Signal');
xlabel('In-Phase'); ylabel('Quadrature');
axis([-2 2 -2 2]); grid on; axis square;

% Plot Transmitted Signal Constellation (reference)
subplot(1, 3, 3);
plot(real(s), imag(s), 'ok');
title('Transmitted 4-QAM Signal');
xlabel('In-Phase'); ylabel('Quadrature');
axis([-2 2 -2 2]); grid on; axis square;

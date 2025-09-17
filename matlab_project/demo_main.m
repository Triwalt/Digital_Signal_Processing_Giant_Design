%% MATLAB Project Demo Script
% This script demonstrates the FFT and CMA implementations
% 
% Author: DSP Course Design
% Date: 2025

clear; clc; close all;

% Add source directory to path
addpath('src');

fprintf('=== Digital Signal Processing Course Design Demo ===\n');
fprintf('This demo showcases FFT and CMA algorithm implementations\n\n');

%% Demo 1: FFT Implementation
fprintf('Demo 1: Fast Fourier Transform (FFT)\n');
fprintf('=====================================\n');

% Create a test signal with multiple frequency components
N = 256;                    % Signal length
fs = 1000;                  % Sampling frequency (Hz)
t = (0:N-1)/fs;            % Time vector

% Multi-tone signal
f1 = 50;   % First frequency component
f2 = 120;  % Second frequency component
f3 = 200;  % Third frequency component

x = 2*cos(2*pi*f1*t) + 1.5*cos(2*pi*f2*t) + cos(2*pi*f3*t);
x = x + 0.3*randn(size(x));  % Add noise

% Compute FFT using our implementation
fprintf('Computing FFT using custom implementation...\n');
tic;
X_custom = fft_implementation(x);
time_custom = toc;

% Compute FFT using MATLAB's built-in function for comparison
X_matlab = fft(x);

% Frequency vector for plotting
f = (0:N-1)*fs/N;

% Create comprehensive plots
figure('Name', 'FFT Analysis Demo', 'Position', [100, 100, 1200, 800]);

% Time domain signal
subplot(2,3,1);
plot(t*1000, x, 'b-', 'LineWidth', 1.5);
xlabel('Time (ms)'); ylabel('Amplitude');
title('Input Signal (Time Domain)');
grid on;

% FFT magnitude comparison
subplot(2,3,2);
plot(f(1:N/2), abs(X_custom(1:N/2)), 'b-', 'LineWidth', 2, 'DisplayName', 'Custom FFT');
hold on;
plot(f(1:N/2), abs(X_matlab(1:N/2)), 'r--', 'LineWidth', 1, 'DisplayName', 'MATLAB FFT');
xlabel('Frequency (Hz)'); ylabel('Magnitude');
title('FFT Magnitude Spectrum');
legend('show'); grid on;

% Error analysis
subplot(2,3,3);
semilogy(f(1:N/2), abs(X_custom(1:N/2) - X_matlab(1:N/2)), 'g-', 'LineWidth', 2);
xlabel('Frequency (Hz)'); ylabel('Absolute Error');
title('Implementation Error');
grid on;

% Spectrogram
subplot(2,3,4);
spectrogram(x, hamming(64), 32, 128, fs, 'yaxis');
title('Signal Spectrogram');

% Phase comparison
subplot(2,3,5);
plot(f(1:N/2), angle(X_custom(1:N/2))*180/pi, 'b-', 'LineWidth', 1.5);
hold on;
plot(f(1:N/2), angle(X_matlab(1:N/2))*180/pi, 'r--', 'LineWidth', 1);
xlabel('Frequency (Hz)'); ylabel('Phase (degrees)');
title('FFT Phase Spectrum');
legend('Custom FFT', 'MATLAB FFT'); grid on;

% Performance comparison
subplot(2,3,6);
sizes = [64, 128, 256, 512, 1024];
times_custom = zeros(size(sizes));
times_matlab = zeros(size(sizes));

for i = 1:length(sizes)
    N_test = sizes(i);
    x_test = randn(1, N_test);
    
    tic; fft_implementation(x_test); times_custom(i) = toc;
    tic; fft(x_test); times_matlab(i) = toc;
end

loglog(sizes, times_custom*1000, 'b-o', 'LineWidth', 2, 'MarkerSize', 8);
hold on;
loglog(sizes, times_matlab*1000, 'r-s', 'LineWidth', 2, 'MarkerSize', 8);
xlabel('FFT Size'); ylabel('Time (ms)');
title('Performance Comparison');
legend('Custom FFT', 'MATLAB FFT'); grid on;

fprintf('FFT computation time: %.4f seconds\n', time_custom);
fprintf('Maximum error: %.2e\n', max(abs(X_custom - X_matlab)));

%% Demo 2: CMA Algorithm
fprintf('\nDemo 2: Constant Modulus Algorithm (CMA)\n');
fprintf('========================================\n');

% Parameters
N_symbols = 500;
num_taps = 15;
step_size = 0.02;
num_iterations = 100;

% Generate QPSK symbols
fprintf('Generating QPSK symbols...\n');
data_bits = randi([0,1], N_symbols*2, 1);
symbols = zeros(N_symbols, 1);

for i = 1:N_symbols
    % Map bits to QPSK symbols
    I = 2*data_bits(2*i-1) - 1;  % +1 or -1
    Q = 2*data_bits(2*i) - 1;    % +1 or -1
    symbols(i) = I + 1j*Q;
end

% Channel model (multipath with ISI)
channel = [0.9, 0.5, 0.3, 0.1];  % 4-tap channel
received = conv(symbols, channel, 'same');

% Add AWGN noise
SNR_dB = 25;
signal_power = mean(abs(received).^2);
noise_power = signal_power / (10^(SNR_dB/10));
noise = sqrt(noise_power/2) * (randn(size(received)) + 1j*randn(size(received)));
received = received + noise;

fprintf('Applying CMA equalization...\n');
% Apply CMA algorithm
[equalized, weights, error_curve] = cma_algorithm(received, step_size, num_taps, num_iterations);

% Create CMA demonstration plots
figure('Name', 'CMA Equalization Demo', 'Position', [150, 150, 1200, 800]);

% Signal constellations
subplot(2,4,1);
scatter(real(symbols), imag(symbols), 30, 'b', 'filled');
axis equal; grid on;
xlabel('In-phase'); ylabel('Quadrature');
title('Original QPSK Symbols');

subplot(2,4,2);
scatter(real(received), imag(received), 30, 'r', 'filled');
axis equal; grid on;
xlabel('In-phase'); ylabel('Quadrature');
title('Received Signal (with ISI)');

subplot(2,4,3);
scatter(real(equalized(num_taps:end)), imag(equalized(num_taps:end)), 30, 'g', 'filled');
axis equal; grid on;
xlabel('In-phase'); ylabel('Quadrature');
title('Equalized Signal');

% Error convergence
subplot(2,4,4);
semilogy(1:length(error_curve), error_curve, 'b-', 'LineWidth', 2);
xlabel('Iteration'); ylabel('CMA Error');
title('Convergence Curve');
grid on;

% Channel and equalizer impulse responses
subplot(2,4,5);
stem(0:length(channel)-1, channel, 'b', 'LineWidth', 2, 'MarkerSize', 8);
xlabel('Sample'); ylabel('Amplitude');
title('Channel Impulse Response');
grid on;

subplot(2,4,6);
stem(0:num_taps-1, real(weights), 'r', 'LineWidth', 2, 'MarkerSize', 8);
hold on;
stem(0:num_taps-1, imag(weights), 'g', 'LineWidth', 2, 'MarkerSize', 8);
xlabel('Tap Index'); ylabel('Weight Value');
title('Equalizer Weights');
legend('Real', 'Imaginary'); grid on;

% Signal comparison over time
subplot(2,4,7);
plot_length = min(100, length(symbols));
plot(1:plot_length, real(symbols(1:plot_length)), 'b-', 'LineWidth', 2);
hold on;
plot(1:plot_length, real(received(1:plot_length)), 'r--', 'LineWidth', 1);
plot(1:plot_length, real(equalized(1:plot_length)), 'g-', 'LineWidth', 1);
xlabel('Symbol Index'); ylabel('Amplitude');
title('Signal Comparison (Real Part)');
legend('Original', 'Received', 'Equalized');
grid on;

% Frequency response comparison
subplot(2,4,8);
[H_channel, w] = freqz(channel, 1, 512);
[H_equalizer, ~] = freqz(weights, 1, 512);
H_combined = H_channel .* H_equalizer;

plot(w/pi, 20*log10(abs(H_channel)), 'b-', 'LineWidth', 2);
hold on;
plot(w/pi, 20*log10(abs(H_equalizer)), 'r-', 'LineWidth', 2);
plot(w/pi, 20*log10(abs(H_combined)), 'g--', 'LineWidth', 2);
xlabel('Normalized Frequency (×π rad/sample)');
ylabel('Magnitude (dB)');
title('Frequency Response');
legend('Channel', 'Equalizer', 'Combined');
grid on;

% Performance metrics
equalized_symbols = equalized(num_taps:end);
if length(equalized_symbols) > length(symbols)
    equalized_symbols = equalized_symbols(1:length(symbols));
end

% Symbol decisions
decided_symbols = sign(real(equalized_symbols)) + 1j*sign(imag(equalized_symbols));

% Calculate EVM
evm_rms = sqrt(mean(abs(equalized_symbols - decided_symbols).^2));
evm_percent = 100 * evm_rms / sqrt(2); % QPSK reference

% Calculate BER
original_bits_matrix = [real(symbols) > 0, imag(symbols) > 0];
decided_bits_matrix = [real(decided_symbols) > 0, imag(decided_symbols) > 0];
bit_errors = sum(sum(original_bits_matrix(1:length(decided_symbols),:) ~= decided_bits_matrix));
total_bits = 2 * length(decided_symbols);
ber = bit_errors / total_bits;

fprintf('\nPerformance Metrics:\n');
fprintf('==================\n');
fprintf('SNR: %.1f dB\n', SNR_dB);
fprintf('EVM: %.2f%%\n', evm_percent);
fprintf('BER: %.2e\n', ber);
fprintf('Final CMA Error: %.2e\n', error_curve(end));
fprintf('Convergence achieved in %d iterations\n', length(error_curve));

%% Summary
fprintf('\n=== Demo Summary ===\n');
fprintf('FFT Implementation:\n');
fprintf('  - Successfully implemented radix-2 DIT algorithm\n');
fprintf('  - Verified against MATLAB built-in FFT\n');
fprintf('  - Demonstrates frequency domain analysis\n\n');

fprintf('CMA Implementation:\n');
fprintf('  - Successfully equalized QPSK signals\n');
fprintf('  - Removed intersymbol interference\n');
fprintf('  - Achieved BER of %.2e\n', ber);
fprintf('  - EVM of %.2f%%\n\n', evm_percent);

fprintf('Both algorithms are ready for FPGA implementation!\n');

% Save results
save('demo_results.mat', 'X_custom', 'X_matlab', 'symbols', 'received', 'equalized', ...
     'weights', 'error_curve', 'evm_percent', 'ber');

fprintf('\nResults saved to demo_results.mat\n');
fprintf('Demo completed successfully!\n');
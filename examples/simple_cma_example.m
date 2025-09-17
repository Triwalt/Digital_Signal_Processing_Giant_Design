%% Simple CMA Example
% Basic example demonstrating CMA equalization
% 
% Author: DSP Course Design
% Date: 2025

clear; clc; close all;

% Add path to source functions
addpath('../matlab_project/src');

fprintf('=== Simple CMA Example ===\n');

% Parameters
N_symbols = 200;           % Number of symbols
num_taps = 11;             % Equalizer taps
step_size = 0.02;          % CMA step size
num_iterations = 30;       % Number of iterations

% Generate QPSK symbols
fprintf('Generating QPSK symbols...\n');
symbols = (2*randi([0,1], N_symbols, 1) - 1) + 1j*(2*randi([0,1], N_symbols, 1) - 1);

% Simple channel model
channel = [1.0, 0.5, 0.2];  % 3-tap channel
received = conv(symbols, channel, 'same');

% Add noise
SNR_dB = 20;
noise_power = mean(abs(received).^2) / (10^(SNR_dB/10));
noise = sqrt(noise_power/2) * (randn(N_symbols,1) + 1j*randn(N_symbols,1));
received = received + noise;

fprintf('Applying CMA equalization...\n');

% Apply CMA
[equalized, weights, error_curve] = cma_algorithm(received, step_size, num_taps, num_iterations);

% Create visualization
figure('Name', 'Simple CMA Example', 'Position', [150, 150, 1000, 600]);

% Signal constellations
subplot(2,3,1);
scatter(real(symbols), imag(symbols), 50, 'b', 'filled');
axis equal; grid on; axis([-2 2 -2 2]);
xlabel('In-phase'); ylabel('Quadrature');
title('Original QPSK');

subplot(2,3,2);
scatter(real(received), imag(received), 50, 'r', 'filled');
axis equal; grid on; axis([-2 2 -2 2]);
xlabel('In-phase'); ylabel('Quadrature');
title('Received (with ISI)');

subplot(2,3,3);
scatter(real(equalized(num_taps:end)), imag(equalized(num_taps:end)), 50, 'g', 'filled');
axis equal; grid on; axis([-2 2 -2 2]);
xlabel('In-phase'); ylabel('Quadrature');
title('Equalized Output');

% Error convergence
subplot(2,3,4);
semilogy(1:length(error_curve), error_curve, 'b-', 'LineWidth', 2);
xlabel('Iteration'); ylabel('CMA Error');
title('Error Convergence');
grid on;

% Equalizer weights
subplot(2,3,5);
stem(1:num_taps, real(weights), 'b', 'LineWidth', 2);
hold on;
stem(1:num_taps, imag(weights), 'r', 'LineWidth', 2);
xlabel('Tap Index'); ylabel('Weight Value');
title('Final Equalizer Weights');
legend('Real', 'Imaginary');
grid on;

% Signal comparison
subplot(2,3,6);
plot_length = min(50, N_symbols);
plot(1:plot_length, real(symbols(1:plot_length)), 'b-', 'LineWidth', 2);
hold on;
plot(1:plot_length, real(received(1:plot_length)), 'r--', 'LineWidth', 1);
plot(1:plot_length, real(equalized(1:plot_length)), 'g-', 'LineWidth', 1);
xlabel('Symbol Index'); ylabel('Amplitude');
title('Signal Comparison');
legend('Original', 'Received', 'Equalized');
grid on;

% Calculate performance metrics
equalized_symbols = equalized(num_taps:end);
decided_symbols = sign(real(equalized_symbols)) + 1j*sign(imag(equalized_symbols));
evm = 100 * sqrt(mean(abs(equalized_symbols - decided_symbols).^2)) / sqrt(2);

fprintf('\nResults:\n');
fprintf('Final CMA Error: %.2e\n', error_curve(end));
fprintf('EVM: %.2f%%\n', evm);
fprintf('Convergence achieved in %d iterations\n', length(error_curve));
fprintf('CMA equalization completed successfully!\n');
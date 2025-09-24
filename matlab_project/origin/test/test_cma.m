%% CMA Algorithm Test Script
% This script tests the Constant Modulus Algorithm implementation
% 
% Author: DSP Course Design
% Date: 2025

clear; clc; close all;

% Add source directory to path
addpath('../src');

%% Test 1: QPSK Signal Equalization
fprintf('=== CMA Algorithm Test ===\n');
fprintf('Test 1: QPSK Signal Equalization\n');

% Parameters
N = 1000;                  % Number of symbols
num_taps = 11;             % Equalizer taps
step_size = 0.01;          % CMA step size
num_iterations = 50;       % Number of CMA iterations

% Generate QPSK symbols
symbols = (2*randi([0,1], N, 1) - 1) + 1j*(2*randi([0,1], N, 1) - 1);

% Channel (introduces ISI)
channel = [0.8, 0.4, 0.2];  % Simple multipath channel
received = conv(symbols, channel);
received = received(1:N);    % Truncate to original length

% Add noise
SNR_dB = 20;
noise_power = mean(abs(received).^2) / (10^(SNR_dB/10));
received = received + sqrt(noise_power/2) * (randn(N,1) + 1j*randn(N,1));

% Apply CMA
fprintf('Running CMA algorithm...\n');
[equalized, weights, error_curve] = cma_algorithm(received, step_size, num_taps, num_iterations);

% Plot results
figure('Name', 'CMA Equalization Results');

% Original and received constellations
subplot(2,3,1);
scatter(real(symbols), imag(symbols), 'b.');
title('Original QPSK Symbols');
xlabel('Real'); ylabel('Imaginary');
axis equal; grid on;

subplot(2,3,2);
scatter(real(received), imag(received), 'r.');
title('Received Signal (with ISI)');
xlabel('Real'); ylabel('Imaginary');
axis equal; grid on;

subplot(2,3,3);
scatter(real(equalized(num_taps:end)), imag(equalized(num_taps:end)), 'g.');
title('Equalized Signal');
xlabel('Real'); ylabel('Imaginary');
axis equal; grid on;

% Error convergence
subplot(2,3,4);
semilogy(1:length(error_curve), error_curve, 'b-', 'LineWidth', 2);
title('CMA Error Convergence');
xlabel('Iteration'); ylabel('Error');
grid on;

% Equalizer weights
subplot(2,3,5);
stem(1:num_taps, real(weights), 'b', 'LineWidth', 2);
hold on;
stem(1:num_taps, imag(weights), 'r', 'LineWidth', 2);
title('Final Equalizer Weights');
xlabel('Tap Index'); ylabel('Weight Value');
legend('Real', 'Imaginary');
grid on;

% Signal comparison over time
subplot(2,3,6);
plot(1:100, real(symbols(1:100)), 'b-', 'LineWidth', 2);
hold on;
plot(1:100, real(received(1:100)), 'r--', 'LineWidth', 1);
plot(1:100, real(equalized(1:100)), 'g-', 'LineWidth', 1);
title('Signal Comparison (Real Part)');
xlabel('Sample'); ylabel('Amplitude');
legend('Original', 'Received', 'Equalized');
grid on;

%% Test 2: Different Step Sizes
fprintf('\nTest 2: Step Size Analysis\n');

step_sizes = [0.001, 0.01, 0.05, 0.1];
colors = {'b', 'r', 'g', 'm'};

figure('Name', 'CMA Step Size Analysis');

for i = 1:length(step_sizes)
    fprintf('Testing step size: %.3f\n', step_sizes(i));
    [~, ~, error_curve_step] = cma_algorithm(received, step_sizes(i), num_taps, num_iterations);
    
    semilogy(1:length(error_curve_step), error_curve_step, ...
             [colors{i} '-'], 'LineWidth', 2, ...
             'DisplayName', sprintf('\\mu = %.3f', step_sizes(i)));
    hold on;
end

title('CMA Convergence for Different Step Sizes');
xlabel('Iteration'); ylabel('Error');
legend('show', 'Location', 'best');
grid on;

%% Test 3: Different Number of Taps
fprintf('\nTest 3: Number of Taps Analysis\n');

tap_counts = [5, 7, 11, 15];
figure('Name', 'CMA Tap Count Analysis');

for i = 1:length(tap_counts)
    fprintf('Testing %d taps\n', tap_counts(i));
    [eq_temp, ~, error_curve_taps] = cma_algorithm(received, 0.01, tap_counts(i), num_iterations);
    
    subplot(2,2,i);
    scatter(real(eq_temp(tap_counts(i):end)), imag(eq_temp(tap_counts(i):end)), 'g.');
    title(sprintf('%d Taps', tap_counts(i)));
    xlabel('Real'); ylabel('Imaginary');
    axis equal; grid on;
end

%% Test 4: Performance Metrics
fprintf('\nTest 4: Performance Metrics\n');

% Calculate EVM (Error Vector Magnitude)
equalized_symbols = equalized(num_taps:end);
% Decision making
decided_symbols = sign(real(equalized_symbols)) + 1j*sign(imag(equalized_symbols));
error_vector = equalized_symbols - decided_symbols;
EVM_percent = 100 * sqrt(mean(abs(error_vector).^2)) / sqrt(mean(abs(decided_symbols).^2));

% Calculate BER (Bit Error Rate) 
original_bits = [real(symbols) > 0, imag(symbols) > 0];
decided_bits = [real(decided_symbols) > 0, imag(decided_symbols) > 0];
bit_errors = sum(sum(original_bits(1:length(decided_symbols),:) ~= decided_bits));
total_bits = 2 * length(decided_symbols);
BER = bit_errors / total_bits;

fprintf('Final Performance Metrics:\n');
fprintf('EVM: %.2f%%\n', EVM_percent);
fprintf('BER: %.2e\n', BER);
fprintf('Final CMA Error: %.2e\n', error_curve(end));
fprintf('Convergence achieved in %d iterations\n', length(error_curve));

%% Test 5: Channel Estimation
fprintf('\nTest 5: Channel Estimation\n');

% Estimate channel response using equalizer weights
estimated_channel = conv(weights, channel);
estimated_channel = estimated_channel / max(abs(estimated_channel));

figure('Name', 'Channel and Equalizer Analysis');

subplot(2,2,1);
stem(1:length(channel), channel, 'b', 'LineWidth', 2);
title('Original Channel');
xlabel('Sample'); ylabel('Amplitude');
grid on;

subplot(2,2,2);
stem(1:num_taps, real(weights), 'r', 'LineWidth', 2);
title('Equalizer Weights (Real)');
xlabel('Tap'); ylabel('Weight');
grid on;

subplot(2,2,3);
f = linspace(0, 1, 512);
H_channel = freqz(channel, 1, f, 'whole');
H_equalizer = freqz(weights, 1, f, 'whole');
plot(f, 20*log10(abs(H_channel)), 'b-', 'LineWidth', 2);
hold on;
plot(f, 20*log10(abs(H_equalizer)), 'r-', 'LineWidth', 2);
plot(f, 20*log10(abs(H_channel .* H_equalizer)), 'g--', 'LineWidth', 2);
title('Frequency Response');
xlabel('Normalized Frequency'); ylabel('Magnitude (dB)');
legend('Channel', 'Equalizer', 'Combined', 'Location', 'best');
grid on;

subplot(2,2,4);
plot(f, angle(H_channel)*180/pi, 'b-', 'LineWidth', 2);
hold on;
plot(f, angle(H_equalizer)*180/pi, 'r-', 'LineWidth', 2);
title('Phase Response');
xlabel('Normalized Frequency'); ylabel('Phase (degrees)');
legend('Channel', 'Equalizer', 'Location', 'best');
grid on;

fprintf('\n=== CMA Test Complete ===\n');
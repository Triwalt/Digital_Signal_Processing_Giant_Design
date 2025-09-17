%% FFT Implementation Test Script
% This script tests the custom FFT implementation against MATLAB's built-in FFT
% 
% Author: DSP Course Design
% Date: 2025

clear; clc; close all;

% Add source directory to path
addpath('../src');

%% Test 1: Simple Sinusoidal Signal
fprintf('=== FFT Implementation Test ===\n');
fprintf('Test 1: Simple Sinusoidal Signal\n');

% Parameters
N = 64;                    % Signal length (power of 2)
fs = 1000;                 % Sampling frequency
t = (0:N-1)/fs;           % Time vector

% Generate test signal: sum of two sinusoids
f1 = 50; f2 = 120;        % Frequencies
x = cos(2*pi*f1*t) + 0.5*cos(2*pi*f2*t);

% Add some noise
x = x + 0.1*randn(size(x));

% Compute FFT using our implementation
tic;
X_custom = fft_implementation(x);
time_custom = toc;

% Compute FFT using MATLAB's built-in function
tic;
X_matlab = fft(x);
time_matlab = toc;

% Compare results
error = max(abs(X_custom - X_matlab));
fprintf('Maximum error between implementations: %.2e\n', error);
fprintf('Custom FFT time: %.4f seconds\n', time_custom);
fprintf('MATLAB FFT time: %.4f seconds\n', time_matlab);

% Plot results
figure('Name', 'FFT Comparison - Sinusoidal Signal');
subplot(2,2,1);
plot(t, real(x));
title('Input Signal (Real Part)');
xlabel('Time (s)'); ylabel('Amplitude');
grid on;

subplot(2,2,2);
f = (0:N-1)*fs/N;
plot(f, abs(X_custom), 'b-', 'LineWidth', 2);
hold on;
plot(f, abs(X_matlab), 'r--', 'LineWidth', 1);
title('FFT Magnitude Comparison');
xlabel('Frequency (Hz)'); ylabel('Magnitude');
legend('Custom FFT', 'MATLAB FFT', 'Location', 'best');
grid on;

subplot(2,2,3);
plot(f, angle(X_custom)*180/pi, 'b-', 'LineWidth', 2);
hold on;
plot(f, angle(X_matlab)*180/pi, 'r--', 'LineWidth', 1);
title('FFT Phase Comparison');
xlabel('Frequency (Hz)'); ylabel('Phase (degrees)');
legend('Custom FFT', 'MATLAB FFT', 'Location', 'best');
grid on;

subplot(2,2,4);
semilogy(f, abs(X_custom - X_matlab));
title('Absolute Error');
xlabel('Frequency (Hz)'); ylabel('Error');
grid on;

%% Test 2: Random Signal
fprintf('\nTest 2: Random Signal\n');

% Generate random signal
x_rand = randn(1, 128) + 1j*randn(1, 128);

% Compute FFTs
X_custom_rand = fft_implementation(x_rand);
X_matlab_rand = fft(x_rand);

% Compare results
error_rand = max(abs(X_custom_rand - X_matlab_rand));
fprintf('Maximum error for random signal: %.2e\n', error_rand);

%% Test 3: Performance Test
fprintf('\nTest 3: Performance Comparison\n');

sizes = [64, 128, 256, 512, 1024];
times_custom = zeros(size(sizes));
times_matlab = zeros(size(sizes));

for i = 1:length(sizes)
    N_test = sizes(i);
    x_test = randn(1, N_test);
    
    % Time custom implementation
    tic;
    for trial = 1:10
        X_test = fft_implementation(x_test);
    end
    times_custom(i) = toc/10;
    
    % Time MATLAB implementation
    tic;
    for trial = 1:10
        X_test_matlab = fft(x_test);
    end
    times_matlab(i) = toc/10;
end

% Plot performance comparison
figure('Name', 'FFT Performance Comparison');
loglog(sizes, times_custom, 'b-o', 'LineWidth', 2, 'MarkerSize', 8);
hold on;
loglog(sizes, times_matlab, 'r-s', 'LineWidth', 2, 'MarkerSize', 8);
xlabel('Signal Length'); ylabel('Time (seconds)');
title('FFT Performance Comparison');
legend('Custom FFT', 'MATLAB FFT', 'Location', 'best');
grid on;

fprintf('Performance test completed.\n');

%% Test 4: Edge Cases
fprintf('\nTest 4: Edge Cases\n');

% Test minimum size (N=2)
x_min = [1, 2];
try
    X_min_custom = fft_implementation(x_min);
    X_min_matlab = fft(x_min);
    error_min = max(abs(X_min_custom - X_min_matlab));
    fprintf('N=2 test passed, error: %.2e\n', error_min);
catch ME
    fprintf('N=2 test failed: %s\n', ME.message);
end

% Test non-power-of-2 (should fail)
x_invalid = [1, 2, 3];
try
    X_invalid = fft_implementation(x_invalid);
    fprintf('Non-power-of-2 test FAILED - should have thrown error\n');
catch ME
    fprintf('Non-power-of-2 test passed - correctly threw error\n');
end

fprintf('\n=== FFT Test Complete ===\n');
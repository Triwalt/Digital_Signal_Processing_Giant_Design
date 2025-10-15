%% Simple FFT Example
% Basic example demonstrating FFT usage
% 
% Author: DSP Course Design
% Date: 2025

clear; clc; close all;

% Add path to source functions
addpath('../matlab_project/Gemini_generated_simulation');

fprintf('=== Simple FFT Example ===\n');

% Parameters
N = 128;                    % Signal length
fs = 1000;                  % Sampling frequency (Hz)
t = (0:N-1)/fs;            % Time vector

% Create a simple test signal
f_signal = 100;             % Signal frequency (Hz)
signal = cos(2*pi*f_signal*t) + 0.3*randn(size(t));

% Compute FFT
fprintf('Computing FFT of %d-point signal...\n', N);
X = my_fft(signal);

% Frequency vector
f = (0:N-1)*fs/N;

% Create plots
figure('Name', 'Simple FFT Example', 'Position', [100, 100, 1000, 600]);

% Time domain
subplot(2,2,1);
plot(t*1000, signal, 'b-', 'LineWidth', 1.5);
xlabel('Time (ms)'); ylabel('Amplitude');
title('Input Signal');
grid on;

% Frequency domain - magnitude
subplot(2,2,2);
plot(f(1:N/2), abs(X(1:N/2)), 'r-', 'LineWidth', 2);
xlabel('Frequency (Hz)'); ylabel('Magnitude');
title('FFT Magnitude Spectrum');
grid on;

% Frequency domain - phase
subplot(2,2,3);
plot(f(1:N/2), angle(X(1:N/2))*180/pi, 'g-', 'LineWidth', 1.5);
xlabel('Frequency (Hz)'); ylabel('Phase (degrees)');
title('FFT Phase Spectrum');
grid on;

% Spectrogram
subplot(2,2,4);
spectrogram(signal, hamming(32), 16, 64, fs, 'yaxis');
title('Signal Spectrogram');

% Find peak frequency
[~, peak_idx] = max(abs(X(1:N/2)));
peak_freq = f(peak_idx);

fprintf('Peak frequency detected at: %.1f Hz (Expected: %.1f Hz)\n', peak_freq, f_signal);
fprintf('FFT computation completed successfully!\n');
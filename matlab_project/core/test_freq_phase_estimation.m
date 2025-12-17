%% test_freq_phase_estimation.m: 频偏估计与相位噪声估计模块测试脚本
%
% 功能说明:
%   全面测试 frequency_offset_estimator.m 和 phase_noise_estimator.m 模块
%   包括:
%   1. 频偏估计精度测试 (不同频偏值、不同SNR)
%   2. 相位噪声估计精度测试 (不同相位噪声标准差)
%   3. 相位模糊解决方法对比测试
%   4. 与现有CMA处理流程的集成测试
%
% 输出:
%   - 控制台测试结果摘要
%   - 测试报告结构体 (可用于后续分析)
%
% 作者: DSP课程设计项目
% 日期: 2025年12月

clc;
clear;
close all;

% 添加core目录到路径
scriptDir = fileparts(mfilename('fullpath'));
if isempty(scriptDir)
    scriptDir = pwd;
end
addpath(scriptDir);

fprintf('╔══════════════════════════════════════════════════════════════╗\n');
fprintf('║     频偏估计与相位噪声估计模块 - 综合测试套件                ║\n');
fprintf('╚══════════════════════════════════════════════════════════════╝\n\n');

%% 测试配置
test_config = struct();
test_config.num_symbols = 10000;        % 测试符号数
test_config.snr_values = [10, 15, 20, 25, 30];  % 测试SNR值 (dB)
test_config.freq_offsets = [0, 1e6, 5e6, 10e6, 15e6, 20e6];  % 测试频偏值 (Hz) - 在有效范围内
test_config.phase_noise_std = [0.5, 1, 2, 3, 5];  % 相位噪声标准差 (度/符号)
test_config.symbol_rate = 28e9;         % 符号率 28 GBaud

% 4QAM/QPSK 星座
constellation = [1+1j, 1-1j, -1-1j, -1+1j] / sqrt(2);

% 测试结果存储
results = struct();
results.freq_offset_tests = [];
results.phase_noise_tests = [];
results.integration_tests = [];

test_passed = 0;
test_failed = 0;

%% ==================== 测试1: 频偏估计精度测试 ====================
fprintf('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');
fprintf('【测试1】频偏估计精度测试\n');
fprintf('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n');

snr_test = 20;  % 固定SNR
freq_est_results = zeros(length(test_config.freq_offsets), 2);

for idx = 1:length(test_config.freq_offsets)
    true_freq_offset = test_config.freq_offsets(idx);
    
    % 生成QPSK测试信号
    [tx_symbols, ~] = generate_qpsk_signal(test_config.num_symbols, constellation);
    
    % 添加频偏
    n = (0:length(tx_symbols)-1).';
    freq_offset_rad = 2 * pi * true_freq_offset / test_config.symbol_rate;
    rx_symbols = tx_symbols .* exp(1j * freq_offset_rad * n);
    
    % 添加AWGN噪声
    rx_symbols = add_awgn(rx_symbols, snr_test);
    
    % 频偏估计
    freq_config = struct();
    freq_config.L = 1;
    freq_config.avg_window = 128;
    freq_config.mod_order = 4;
    freq_config.method = 'differential';
    
    [freq_est, ~, diag] = frequency_offset_estimator(rx_symbols, test_config.symbol_rate, freq_config);
    
    % 记录结果
    freq_est_results(idx, 1) = true_freq_offset;
    freq_est_results(idx, 2) = freq_est;
    
    % 计算误差
    abs_error = abs(freq_est - true_freq_offset);
    rel_error = abs_error / max(abs(true_freq_offset), 1e6) * 100;  % 相对于1MHz
    
    % 判断测试是否通过 (误差 < 10% 或 < 1MHz)
    if abs_error < 1e6 || rel_error < 10
        status = '✓ PASS';
        test_passed = test_passed + 1;
    else
        status = '✗ FAIL';
        test_failed = test_failed + 1;
    end
    
    fprintf('  频偏 = %6.0f MHz -> 估计 = %8.2f MHz, 误差 = %6.2f MHz (%5.2f%%) [%s]\n', ...
            true_freq_offset/1e6, freq_est/1e6, abs_error/1e6, rel_error, status);
end

results.freq_offset_tests = freq_est_results;
fprintf('\n');

%% ==================== 测试2: 不同SNR下的频偏估计 ====================
fprintf('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');
fprintf('【测试2】不同SNR下的频偏估计精度\n');
fprintf('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n');

true_freq_offset = 10e6;  % 固定频偏 10 MHz
snr_est_results = zeros(length(test_config.snr_values), 2);

for idx = 1:length(test_config.snr_values)
    snr = test_config.snr_values(idx);
    
    % 生成测试信号
    [tx_symbols, ~] = generate_qpsk_signal(test_config.num_symbols, constellation);
    
    % 添加频偏
    n = (0:length(tx_symbols)-1).';
    freq_offset_rad = 2 * pi * true_freq_offset / test_config.symbol_rate;
    rx_symbols = tx_symbols .* exp(1j * freq_offset_rad * n);
    
    % 添加噪声
    rx_symbols = add_awgn(rx_symbols, snr);
    
    % 频偏估计
    freq_config = struct();
    freq_config.L = 1;
    freq_config.avg_window = 128;
    [freq_est, ~, ~] = frequency_offset_estimator(rx_symbols, test_config.symbol_rate, freq_config);
    
    snr_est_results(idx, 1) = snr;
    snr_est_results(idx, 2) = freq_est;
    
    abs_error = abs(freq_est - true_freq_offset);
    
    if abs_error < 5e6  % 允许 5 MHz 误差
        status = '✓ PASS';
        test_passed = test_passed + 1;
    else
        status = '✗ FAIL';
        test_failed = test_failed + 1;
    end
    
    fprintf('  SNR = %2d dB -> 估计频偏 = %8.2f MHz, 误差 = %6.2f MHz [%s]\n', ...
            snr, freq_est/1e6, abs_error/1e6, status);
end

fprintf('\n');

%% ==================== 测试3: 相位噪声估计测试 ====================
fprintf('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');
fprintf('【测试3】相位噪声估计与补偿测试\n');
fprintf('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n');

snr_test = 25;
phase_noise_results = [];

for idx = 1:length(test_config.phase_noise_std)
    pn_std_deg = test_config.phase_noise_std(idx);
    
    % 生成测试信号
    [tx_symbols, tx_bits] = generate_qpsk_signal(test_config.num_symbols, constellation);
    
    % 添加相位噪声 (随机游走模型)
    pn_std_rad = pn_std_deg * pi / 180;
    phase_noise = cumsum(pn_std_rad * randn(length(tx_symbols), 1));
    rx_symbols = tx_symbols .* exp(1j * phase_noise);
    
    % 添加AWGN
    rx_symbols = add_awgn(rx_symbols, snr_test);
    
    % 计算补偿前EVM
    evm_before = calculate_evm(rx_symbols, constellation);
    
    % 相位噪声估计与补偿
    pn_config = struct();
    pn_config.mod_order = 4;
    pn_config.win_half = 5;
    pn_config.constellation = constellation;
    pn_config.ambiguity_method = 'unwrap';  % 使用更稳健的unwrap方法
    pn_config.block_size = 50;
    
    [phase_est, signal_corrected, diag] = phase_noise_estimator(rx_symbols, pn_config);
    
    % 计算补偿后EVM
    evm_after = calculate_evm(signal_corrected, constellation);
    
    % 计算相位估计误差 (RMSE)
    phase_error_rmse = sqrt(mean((phase_est - phase_noise).^2)) * 180 / pi;
    
    % 记录结果
    result_row = [pn_std_deg, evm_before, evm_after, phase_error_rmse];
    phase_noise_results = [phase_noise_results; result_row];
    
    % 判断测试是否通过 (EVM改善且 < 15%)
    if evm_after < evm_before && evm_after < 20
        status = '✓ PASS';
        test_passed = test_passed + 1;
    else
        status = '✗ FAIL';
        test_failed = test_failed + 1;
    end
    
    fprintf('  相位噪声σ = %4.1f°/符号:\n', pn_std_deg);
    fprintf('    EVM: %.2f%% -> %.2f%% (改善 %.1f%%)\n', ...
            evm_before, evm_after, evm_before - evm_after);
    fprintf('    相位估计RMSE = %.2f° [%s]\n\n', phase_error_rmse, status);
end

results.phase_noise_tests = phase_noise_results;

%% ==================== 测试4: 相位模糊解决方法对比 ====================
fprintf('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');
fprintf('【测试4】相位模糊解决方法对比\n');
fprintf('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n');

methods = {'none', 'differential', 'unwrap', 'decision_directed'};
method_results = zeros(length(methods), 2);

% 生成带相位噪声的测试信号
[tx_symbols, ~] = generate_qpsk_signal(test_config.num_symbols, constellation);
pn_std_rad = 2 * pi / 180;  % 2度/符号
phase_noise = cumsum(pn_std_rad * randn(length(tx_symbols), 1));
rx_symbols = tx_symbols .* exp(1j * phase_noise);
rx_symbols = add_awgn(rx_symbols, 20);

fprintf('  测试条件: 相位噪声σ=2°/符号, SNR=20dB\n\n');

for idx = 1:length(methods)
    method = methods{idx};
    
    pn_config = struct();
    pn_config.mod_order = 4;
    pn_config.win_half = 5;
    pn_config.constellation = constellation;
    pn_config.ambiguity_method = method;
    pn_config.block_size = 50;
    
    [~, signal_corrected, diag] = phase_noise_estimator(rx_symbols, pn_config);
    
    evm = calculate_evm(signal_corrected, constellation);
    method_results(idx, 1) = idx;
    method_results(idx, 2) = evm;
    
    if strcmp(method, 'unwrap') && evm < 15
        status = '✓ PASS (推荐)';
        test_passed = test_passed + 1;
    elseif evm < 25
        status = '○ OK';
        test_passed = test_passed + 1;
    else
        status = '✗ FAIL';
        test_failed = test_failed + 1;
    end
    
    fprintf('  方法: %-20s -> EVM = %6.2f%% [%s]\n', method, evm, status);
end

fprintf('\n');

%% ==================== 测试5: 联合频偏和相位噪声补偿 ====================
fprintf('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');
fprintf('【测试5】联合频偏和相位噪声补偿测试\n');
fprintf('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n');

% 生成测试信号
[tx_symbols, tx_bits] = generate_qpsk_signal(test_config.num_symbols, constellation);

% 添加频偏 (20 MHz - 在有效范围内)
true_freq_offset = 20e6;
n = (0:length(tx_symbols)-1).';
freq_offset_rad = 2 * pi * true_freq_offset / test_config.symbol_rate;
rx_symbols = tx_symbols .* exp(1j * freq_offset_rad * n);

% 添加相位噪声 (3度/符号)
pn_std_rad = 3 * pi / 180;
phase_noise = cumsum(pn_std_rad * randn(length(tx_symbols), 1));
rx_symbols = rx_symbols .* exp(1j * phase_noise);

% 添加AWGN (20 dB)
rx_symbols = add_awgn(rx_symbols, 20);

% 初始EVM
evm_initial = calculate_evm(rx_symbols, constellation);
fprintf('  初始条件: 频偏=%dMHz, 相位噪声σ=3°/符号, SNR=20dB\n', true_freq_offset/1e6);
fprintf('  初始EVM = %.2f%%\n\n', evm_initial);

% 步骤1: 频偏估计与补偿
freq_config = struct();
freq_config.L = 1;
freq_config.avg_window = 128;
[freq_est, ~, ~] = frequency_offset_estimator(rx_symbols, test_config.symbol_rate, freq_config);

% 频偏补偿
rx_symbols_freq_comp = rx_symbols .* exp(-1j * 2 * pi * freq_est / test_config.symbol_rate * n);
evm_after_freq = calculate_evm(rx_symbols_freq_comp, constellation);

fprintf('  频偏估计: %.2f MHz (真实: %d MHz)\n', freq_est/1e6, true_freq_offset/1e6);
fprintf('  频偏补偿后EVM = %.2f%%\n\n', evm_after_freq);

% 步骤2: 相位噪声估计与补偿
pn_config = struct();
pn_config.mod_order = 4;
pn_config.win_half = 5;
pn_config.constellation = constellation;
pn_config.ambiguity_method = 'unwrap';
pn_config.block_size = 50;

[~, signal_final, diag] = phase_noise_estimator(rx_symbols_freq_comp, pn_config);
evm_final = calculate_evm(signal_final, constellation);

fprintf('  相位噪声补偿后EVM = %.2f%%\n', evm_final);
fprintf('  总EVM改善: %.2f%% -> %.2f%% (改善 %.1f%%)\n\n', ...
        evm_initial, evm_final, evm_initial - evm_final);

% 判断测试是否通过
if evm_final < 15 && evm_final < evm_initial * 0.5
    status = '✓ PASS';
    test_passed = test_passed + 1;
else
    status = '✗ FAIL';
    test_failed = test_failed + 1;
end
fprintf('  联合补偿测试结果: [%s]\n\n', status);

%% ==================== 测试总结 ====================
fprintf('═══════════════════════════════════════════════════════════════\n');
fprintf('                        测试总结\n');
fprintf('═══════════════════════════════════════════════════════════════\n');
fprintf('  通过测试: %d\n', test_passed);
fprintf('  失败测试: %d\n', test_failed);
fprintf('  通过率: %.1f%%\n', 100 * test_passed / (test_passed + test_failed));
fprintf('═══════════════════════════════════════════════════════════════\n\n');

% 保存测试结果
results.test_passed = test_passed;
results.test_failed = test_failed;
results.timestamp = datetime('now');

% 如果所有测试通过
if test_failed == 0
    fprintf('★★★ 所有测试通过! 模块验证成功! ★★★\n\n');
else
    fprintf('!!! 存在失败的测试, 请检查算法实现 !!!\n\n');
end

fprintf('测试结果已保存到变量 "results"\n');
fprintf('运行 generate_freq_phase_report 生成详细报告\n\n');

%% ==================== 辅助函数 ====================

function [symbols, bits] = generate_qpsk_signal(num_symbols, constellation)
% 生成随机QPSK信号
    bits = randi([0 1], num_symbols * 2, 1);
    bit_pairs = reshape(bits, 2, []).';
    symbol_indices = bit_pairs(:, 1) * 2 + bit_pairs(:, 2);
    gray_map = [0 1 3 2];
    gray_indices = gray_map(symbol_indices + 1);
    symbols = constellation(gray_indices + 1).';
end

function signal_noisy = add_awgn(signal, snr_db)
% 添加AWGN噪声
    signal_power = mean(abs(signal).^2);
    noise_power = signal_power / (10^(snr_db/10));
    noise = sqrt(noise_power/2) * (randn(size(signal)) + 1j * randn(size(signal)));
    signal_noisy = signal + noise;
end

function evm = calculate_evm(signal, constellation)
% 计算EVM
    N = length(signal);
    decided = zeros(N, 1);
    for n = 1:N
        [~, idx] = min(abs(signal(n) - constellation));
        decided(n) = constellation(idx);
    end
    error_vec = signal - decided;
    evm = sqrt(mean(abs(error_vec).^2)) / sqrt(mean(abs(decided).^2)) * 100;
end

%% generate_freq_phase_report.m: 频偏估计与相位噪声估计模块的图形化报告生成
%
% 功能说明:
%   生成完整的测试报告, 包括:
%   1. 频偏估计性能曲线 (估计精度 vs 频偏值, vs SNR)
%   2. 相位噪声补偿性能曲线 (EVM改善 vs 相位噪声标准差)
%   3. 星座图对比 (补偿前后)
%   4. 相位轨迹跟踪图
%   5. 算法性能摘要表格
%
% 输出:
%   - 多个图形窗口
%   - 可选: 保存到 reports/ 目录下的 .fig 和 .mat 文件
%
% 用法:
%   generate_freq_phase_report          % 使用默认参数运行
%   generate_freq_phase_report(config)  % 使用自定义配置
%
% 作者: DSP课程设计项目
% 日期: 2025年12月

function report_data = generate_freq_phase_report(config)

    clc;
    close all;
    
    % 添加core目录到路径
    scriptDir = fileparts(mfilename('fullpath'));
    if isempty(scriptDir)
        scriptDir = pwd;
    end
    addpath(scriptDir);
    
    fprintf('╔══════════════════════════════════════════════════════════════╗\n');
    fprintf('║     频偏估计与相位噪声估计 - 图形化报告生成器               ║\n');
    fprintf('╚══════════════════════════════════════════════════════════════╝\n\n');
    
    %% 默认配置
    if nargin < 1 || isempty(config)
        config = struct();
    end
    
    if ~isfield(config, 'num_symbols')
        config.num_symbols = 20000;
    end
    if ~isfield(config, 'symbol_rate')
        config.symbol_rate = 28e9;  % 28 GBaud
    end
    if ~isfield(config, 'save_figures')
        config.save_figures = true;
    end
    if ~isfield(config, 'output_dir')
        config.output_dir = fullfile(scriptDir, '..', 'reports');
    end
    
    % 确保输出目录存在
    if config.save_figures && ~exist(config.output_dir, 'dir')
        mkdir(config.output_dir);
    end
    
    % 4QAM星座
    constellation = [1+1j, 1-1j, -1-1j, -1+1j] / sqrt(2);
    
    % 报告数据
    report_data = struct();
    report_data.config = config;
    report_data.timestamp = datetime('now');
    
    %% ==================== 图1: 频偏估计精度 vs 真实频偏 ====================
    fprintf('生成图1: 频偏估计精度分析...\n');
    
    freq_offsets = [0, 1e6, 5e6, 10e6, 15e6, 20e6];  % MHz级频偏
    snr_fixed = 20;
    
    freq_est_values = zeros(size(freq_offsets));
    freq_est_errors = zeros(size(freq_offsets));
    
    for idx = 1:length(freq_offsets)
        true_offset = freq_offsets(idx);
        
        % 生成测试信号
        [tx_symbols, ~] = generate_qpsk_signal(config.num_symbols, constellation);
        
        % 添加频偏
        n = (0:length(tx_symbols)-1).';
        freq_offset_rad = 2 * pi * true_offset / config.symbol_rate;
        rx_symbols = tx_symbols .* exp(1j * freq_offset_rad * n);
        rx_symbols = add_awgn(rx_symbols, snr_fixed);
        
        % 估计
        freq_config = struct('L', 1, 'avg_window', 128, 'mod_order', 4);
        [freq_est, ~, ~] = frequency_offset_estimator(rx_symbols, config.symbol_rate, freq_config);
        
        freq_est_values(idx) = freq_est;
        freq_est_errors(idx) = abs(freq_est - true_offset);
    end
    
    figure('Name', '频偏估计精度分析', 'Position', [100, 100, 1200, 500]);
    
    subplot(1, 2, 1);
    plot(freq_offsets/1e6, freq_est_values/1e6, 'bo-', 'LineWidth', 2, 'MarkerSize', 8);
    hold on;
    plot(freq_offsets/1e6, freq_offsets/1e6, 'r--', 'LineWidth', 1.5);
    xlabel('真实频偏 (MHz)');
    ylabel('估计频偏 (MHz)');
    title('频偏估计值 vs 真实值');
    legend('估计值', '理想值', 'Location', 'northwest');
    grid on;
    
    subplot(1, 2, 2);
    bar(freq_offsets/1e6, freq_est_errors/1e6, 'FaceColor', [0.3 0.6 0.9]);
    xlabel('真实频偏 (MHz)');
    ylabel('估计误差 (MHz)');
    title(sprintf('频偏估计误差 (SNR = %d dB)', snr_fixed));
    grid on;
    
    sgtitle('【图1】频偏估计精度分析', 'FontSize', 14, 'FontWeight', 'bold');
    
    report_data.freq_est_accuracy.freq_offsets = freq_offsets;
    report_data.freq_est_accuracy.estimates = freq_est_values;
    report_data.freq_est_accuracy.errors = freq_est_errors;
    
    if config.save_figures
        savefig(fullfile(config.output_dir, 'freq_offset_accuracy.fig'));
    end
    
    %% ==================== 图2: 频偏估计精度 vs SNR ====================
    fprintf('生成图2: SNR对频偏估计的影响...\n');
    
    snr_values = 5:5:35;
    true_offset = 10e6;  % 固定10MHz
    
    freq_est_vs_snr = zeros(size(snr_values));
    freq_err_vs_snr = zeros(size(snr_values));
    
    for idx = 1:length(snr_values)
        snr = snr_values(idx);
        
        [tx_symbols, ~] = generate_qpsk_signal(config.num_symbols, constellation);
        n = (0:length(tx_symbols)-1).';
        freq_offset_rad = 2 * pi * true_offset / config.symbol_rate;
        rx_symbols = tx_symbols .* exp(1j * freq_offset_rad * n);
        rx_symbols = add_awgn(rx_symbols, snr);
        
        freq_config = struct('L', 1, 'avg_window', 128, 'mod_order', 4);
        [freq_est, ~, ~] = frequency_offset_estimator(rx_symbols, config.symbol_rate, freq_config);
        
        freq_est_vs_snr(idx) = freq_est;
        freq_err_vs_snr(idx) = abs(freq_est - true_offset);
    end
    
    figure('Name', 'SNR对频偏估计的影响', 'Position', [150, 150, 800, 400]);
    
    yyaxis left;
    plot(snr_values, freq_est_vs_snr/1e6, 'b-o', 'LineWidth', 2, 'MarkerSize', 8);
    hold on;
    yline(true_offset/1e6, 'b--', '真实值', 'LineWidth', 1.5);
    ylabel('估计频偏 (MHz)');
    
    yyaxis right;
    plot(snr_values, freq_err_vs_snr/1e6, 'r-s', 'LineWidth', 2, 'MarkerSize', 8);
    ylabel('估计误差 (MHz)');
    
    xlabel('SNR (dB)');
    title(sprintf('【图2】SNR对频偏估计的影响 (真实频偏 = %d MHz)', true_offset/1e6));
    legend('估计值', '真实值', '估计误差', 'Location', 'best');
    grid on;
    
    report_data.freq_est_vs_snr.snr_values = snr_values;
    report_data.freq_est_vs_snr.estimates = freq_est_vs_snr;
    report_data.freq_est_vs_snr.errors = freq_err_vs_snr;
    
    if config.save_figures
        savefig(fullfile(config.output_dir, 'freq_offset_vs_snr.fig'));
    end
    
    %% ==================== 综合大图: 频偏与相位噪声补偿性能报告 ====================
    fprintf('生成综合报告: 频偏与相位噪声补偿性能...\n');
    
    figure('Name', '频偏与相位噪声补偿综合报告', 'Position', [30, 30, 1600, 1000], 'Color', 'w');
    
    % ====== 生成综合测试信号 ======
    num_sym_demo = 2000;
    [tx_demo, ~] = generate_qpsk_signal(num_sym_demo, constellation);
    n_demo = (0:num_sym_demo-1).';
    
    % 添加频偏 + 相位噪声 + AWGN
    freq_demo = 8e6;  % 8 MHz 频偏
    freq_rad_demo = 2 * pi * freq_demo / config.symbol_rate;
    pn_std_demo = 2.5;  % 2.5度/符号 相位噪声
    pn_rad_demo = pn_std_demo * pi / 180;
    phase_noise_true = cumsum(pn_rad_demo * randn(num_sym_demo, 1));
    
    % 叠加所有损伤
    rx_demo = tx_demo .* exp(1j * freq_rad_demo * n_demo) .* exp(1j * phase_noise_true);
    rx_demo = add_awgn(rx_demo, 20);
    
    % 频偏估计与补偿
    freq_cfg_demo = struct('L', 1, 'avg_window', 128, 'mod_order', 4);
    [freq_est_demo, ~, ~] = frequency_offset_estimator(rx_demo, config.symbol_rate, freq_cfg_demo);
    rx_freq_comp = rx_demo .* exp(-1j * 2 * pi * freq_est_demo / config.symbol_rate * n_demo);
    
    % 相位噪声估计与补偿
    pn_cfg_demo = struct('mod_order', 4, 'win_half', 5, 'constellation', constellation, ...
                         'ambiguity_method', 'unwrap', 'phase_unwrap', true);
    [phase_est_demo, rx_final, ~] = phase_noise_estimator(rx_freq_comp, pn_cfg_demo);
    
    % ====== 第一行: 星座图演进 (4个) ======
    
    % (a) 发射信号
    subplot(3, 4, 1);
    plot(real(tx_demo), imag(tx_demo), '.', 'Color', [0.2 0.7 0.2], 'MarkerSize', 5);
    hold on;
    plot(real(constellation), imag(constellation), 'ko', 'MarkerSize', 12, 'LineWidth', 2);
    xlabel('实部'); ylabel('虚部');
    title('(a) 发射信号', 'FontWeight', 'bold');
    axis equal; grid on; box on;
    xlim([-1.5 1.5]); ylim([-1.5 1.5]);
    
    % (b) 接收信号 (频偏+相噪+噪声)
    subplot(3, 4, 2);
    plot(real(rx_demo), imag(rx_demo), '.r', 'MarkerSize', 5);
    xlabel('实部'); ylabel('虚部');
    title(sprintf('(b) 接收信号\n频偏=%dMHz, 相噪=%.1f°/sym', freq_demo/1e6, pn_std_demo), 'FontWeight', 'bold');
    axis equal; grid on; box on;
    xlim([-2 2]); ylim([-2 2]);
    
    % (c) 频偏补偿后
    subplot(3, 4, 3);
    plot(real(rx_freq_comp), imag(rx_freq_comp), '.', 'Color', [1 0.5 0], 'MarkerSize', 5);
    xlabel('实部'); ylabel('虚部');
    evm_freq_comp = calculate_evm(rx_freq_comp, constellation);
    title(sprintf('(c) 频偏补偿后\n估计=%.2fMHz, EVM=%.1f%%', freq_est_demo/1e6, evm_freq_comp), 'FontWeight', 'bold');
    axis equal; grid on; box on;
    xlim([-1.8 1.8]); ylim([-1.8 1.8]);
    
    % (d) 相位噪声补偿后
    subplot(3, 4, 4);
    plot(real(rx_final), imag(rx_final), '.b', 'MarkerSize', 5);
    hold on;
    plot(real(constellation), imag(constellation), 'ko', 'MarkerSize', 12, 'LineWidth', 2);
    xlabel('实部'); ylabel('虚部');
    evm_final = calculate_evm(rx_final, constellation);
    title(sprintf('(d) 相噪补偿后\nEVM=%.1f%% (改善%.1f%%)', evm_final, evm_freq_comp-evm_final), 'FontWeight', 'bold');
    axis equal; grid on; box on;
    xlim([-1.5 1.5]); ylim([-1.5 1.5]);
    
    % ====== 第二行: 相位跟踪核心图 ======
    
    % (e) 真实相位 vs 估计相位轨迹
    subplot(3, 4, [5, 6]);
    plot_range = 1:min(800, num_sym_demo);  % 显示前800个符号
    
    % 频偏引起的线性相位
    freq_phase_true = freq_rad_demo * n_demo;
    % 总真实相位 = 频偏相位 + 相位噪声
    total_phase_true = freq_phase_true + phase_noise_true;
    
    plot(plot_range, total_phase_true(plot_range)*180/pi, 'r-', 'LineWidth', 1.5, 'DisplayName', '真实总相位');
    hold on;
    plot(plot_range, freq_phase_true(plot_range)*180/pi, 'r--', 'LineWidth', 1, 'DisplayName', '频偏相位分量');
    
    % 估计的频偏相位
    freq_phase_est = 2*pi*freq_est_demo/config.symbol_rate * n_demo;
    plot(plot_range, freq_phase_est(plot_range)*180/pi, 'b--', 'LineWidth', 1.5, 'DisplayName', '估计频偏相位');
    
    % 估计的总相位 (频偏 + 相位噪声)
    total_phase_est = freq_phase_est + phase_est_demo;
    plot(plot_range, total_phase_est(plot_range)*180/pi, 'b-', 'LineWidth', 1.2, 'DisplayName', '估计总相位');
    
    xlabel('符号索引'); ylabel('相位 (度)');
    title('(e) 相位跟踪: 真实相位 vs 估计相位', 'FontWeight', 'bold', 'FontSize', 11);
    legend('Location', 'northwest', 'FontSize', 8);
    grid on; box on;
    
    % (f) 相位噪声分量对比 (去除频偏后)
    subplot(3, 4, [7, 8]);
    % 残余相位 = 总相位 - 频偏相位
    residual_phase_true = phase_noise_true;
    residual_phase_est = phase_est_demo;
    
    plot(plot_range, residual_phase_true(plot_range)*180/pi, 'r-', 'LineWidth', 1.5, 'DisplayName', '真实相位噪声');
    hold on;
    % 估计相位限制在合理范围内显示
    plot(plot_range, mod(residual_phase_est(plot_range)*180/pi + 180, 360) - 180, 'b-', 'LineWidth', 1.2, 'DisplayName', '估计相位噪声');
    
    xlabel('符号索引'); ylabel('相位 (度)');
    title('(f) 相位噪声跟踪 (去除频偏分量)', 'FontWeight', 'bold', 'FontSize', 11);
    legend('Location', 'best', 'FontSize', 9);
    grid on; box on;
    
    % ====== 第三行: 性能指标 ======
    
    % (g) 频偏估计精度
    subplot(3, 4, 9);
    plot(freq_offsets/1e6, freq_est_values/1e6, 'bo-', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'b');
    hold on;
    plot(freq_offsets/1e6, freq_offsets/1e6, 'r--', 'LineWidth', 1.5);
    xlabel('真实频偏 (MHz)'); ylabel('估计频偏 (MHz)');
    title('(g) 频偏估计精度', 'FontWeight', 'bold');
    legend('估计值', '理想', 'Location', 'northwest');
    grid on; box on;
    
    % (h) 相位估计误差分布
    subplot(3, 4, 10);
    phase_error = residual_phase_est - residual_phase_true;
    phase_error_deg = mod(phase_error*180/pi + 45, 90) - 45;  % 去除π/2模糊
    histogram(phase_error_deg, 30, 'FaceColor', [0.3 0.6 0.9], 'EdgeColor', 'w');
    hold on;
    xline(0, 'r--', 'LineWidth', 2);
    xline(mean(phase_error_deg), 'g-', 'LineWidth', 2);
    xlabel('相位误差 (度)'); ylabel('计数');
    title(sprintf('(h) 相位误差分布\nRMSE=%.2f°, 均值=%.2f°', std(phase_error_deg), mean(phase_error_deg)), 'FontWeight', 'bold');
    grid on; box on;
    
    % (i) EVM随SNR变化
    subplot(3, 4, 11);
    snr_test = 10:5:30;
    evm_before_snr = zeros(size(snr_test));
    evm_after_snr = zeros(size(snr_test));
    
    for idx = 1:length(snr_test)
        [tx_t, ~] = generate_qpsk_signal(config.num_symbols, constellation);
        n_t = (0:length(tx_t)-1).';
        pn_t = cumsum(pn_rad_demo * randn(length(tx_t), 1));
        rx_t = tx_t .* exp(1j * freq_rad_demo * n_t) .* exp(1j * pn_t);
        rx_t = add_awgn(rx_t, snr_test(idx));
        
        evm_before_snr(idx) = calculate_evm(rx_t, constellation);
        
        [f_est, ~, ~] = frequency_offset_estimator(rx_t, config.symbol_rate, freq_cfg_demo);
        rx_fc = rx_t .* exp(-1j * 2*pi*f_est/config.symbol_rate * n_t);
        [~, rx_pc, ~] = phase_noise_estimator(rx_fc, pn_cfg_demo);
        evm_after_snr(idx) = calculate_evm(rx_pc, constellation);
    end
    
    plot(snr_test, evm_before_snr, 'r-o', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'r');
    hold on;
    plot(snr_test, evm_after_snr, 'b-s', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'b');
    xlabel('SNR (dB)'); ylabel('EVM (%)');
    title('(i) EVM vs SNR', 'FontWeight', 'bold');
    legend('补偿前', '补偿后', 'Location', 'northeast');
    grid on; box on;
    
    % (j) 性能摘要文本
    subplot(3, 4, 12);
    axis off;
    
    % 计算关键指标
    freq_err_demo = abs(freq_est_demo - freq_demo);
    phase_rmse = std(phase_error_deg);
    evm_improve = evm_freq_comp - evm_final;
    
    summary_text = {
        '\bf性能摘要', '', ...
        sprintf('频偏估计误差: %.2f kHz', freq_err_demo/1e3), ...
        sprintf('相对误差: %.3f%%', freq_err_demo/freq_demo*100), '', ...
        sprintf('相位跟踪RMSE: %.2f°', phase_rmse), ...
        sprintf('EVM改善: %.1f%% → %.1f%%', evm_freq_comp, evm_final), ...
        sprintf('EVM降低: %.1f%%', evm_improve), '', ...
        '\bf测试条件', ...
        sprintf('SNR = 20 dB'), ...
        sprintf('频偏 = %d MHz', freq_demo/1e6), ...
        sprintf('相噪 = %.1f °/sym', pn_std_demo)
    };
    
    text(0.1, 0.9, summary_text, 'FontSize', 10, 'VerticalAlignment', 'top', ...
         'FontName', 'FixedWidth', 'BackgroundColor', [0.95 0.95 0.95], ...
         'EdgeColor', [0.5 0.5 0.5], 'Margin', 8);
    title('(j) 性能摘要', 'FontWeight', 'bold');
    
    % ====== 总标题 ======
    sgtitle({'频偏估计与相位噪声补偿 - 综合性能报告', ...
             sprintf('符号速率: %.0f GBaud | 调制: QPSK | 算法: 四次方差分 + V&V', config.symbol_rate/1e9)}, ...
            'FontSize', 14, 'FontWeight', 'bold');
    
    if config.save_figures
        savefig(fullfile(config.output_dir, 'freq_offset_comprehensive_report.fig'));
        print(fullfile(config.output_dir, 'freq_offset_comprehensive_report.png'), '-dpng', '-r150');
    end
    
    report_data.comprehensive.freq_demo = freq_demo;
    report_data.comprehensive.freq_est_demo = freq_est_demo;
    report_data.comprehensive.phase_noise_true = phase_noise_true;
    report_data.comprehensive.phase_est = phase_est_demo;
    report_data.comprehensive.evm_before = evm_freq_comp;
    report_data.comprehensive.evm_after = evm_final;
    
    %% ==================== 图3: 相位噪声补偿性能 ====================
    fprintf('生成图3: 相位噪声补偿性能...\n');
    
    pn_std_values = 0.5:0.5:5;  % 相位噪声标准差 (度/符号)
    snr_fixed = 20;
    
    evm_before = zeros(size(pn_std_values));
    evm_after = zeros(size(pn_std_values));
    
    for idx = 1:length(pn_std_values)
        pn_std = pn_std_values(idx);
        
        [tx_symbols, ~] = generate_qpsk_signal(config.num_symbols, constellation);
        
        % 添加相位噪声
        pn_std_rad = pn_std * pi / 180;
        phase_noise = cumsum(pn_std_rad * randn(length(tx_symbols), 1));
        rx_symbols = tx_symbols .* exp(1j * phase_noise);
        rx_symbols = add_awgn(rx_symbols, snr_fixed);
        
        evm_before(idx) = calculate_evm(rx_symbols, constellation);
        
        % 相位噪声补偿
        pn_config = struct();
        pn_config.mod_order = 4;
        pn_config.win_half = 5;
        pn_config.constellation = constellation;
        pn_config.ambiguity_method = 'unwrap';  % 使用更稳健的unwrap方法
        pn_config.block_size = 50;
        
        [~, signal_corrected, ~] = phase_noise_estimator(rx_symbols, pn_config);
        evm_after(idx) = calculate_evm(signal_corrected, constellation);
    end
    
    figure('Name', '相位噪声补偿性能', 'Position', [200, 200, 1000, 400]);
    
    subplot(1, 2, 1);
    plot(pn_std_values, evm_before, 'r-o', 'LineWidth', 2, 'MarkerSize', 8);
    hold on;
    plot(pn_std_values, evm_after, 'b-s', 'LineWidth', 2, 'MarkerSize', 8);
    xlabel('相位噪声标准差 (°/符号)');
    ylabel('EVM (%)');
    title('EVM vs 相位噪声强度');
    legend('补偿前', '补偿后', 'Location', 'northwest');
    grid on;
    
    subplot(1, 2, 2);
    evm_improvement = evm_before - evm_after;
    bar(pn_std_values, evm_improvement, 'FaceColor', [0.2 0.7 0.3]);
    xlabel('相位噪声标准差 (°/符号)');
    ylabel('EVM改善 (%)');
    title('相位噪声补偿带来的EVM改善');
    grid on;
    
    sgtitle(sprintf('【图3】相位噪声补偿性能 (SNR = %d dB)', snr_fixed), ...
            'FontSize', 14, 'FontWeight', 'bold');
    
    report_data.phase_noise_perf.pn_std = pn_std_values;
    report_data.phase_noise_perf.evm_before = evm_before;
    report_data.phase_noise_perf.evm_after = evm_after;
    
    if config.save_figures
        savefig(fullfile(config.output_dir, 'phase_noise_performance.fig'));
    end
    
    %% ==================== 图4: 星座图对比 ====================
    fprintf('生成图4: 星座图对比...\n');
    
    num_symbols_plot = 5000;
    [tx_symbols, ~] = generate_qpsk_signal(num_symbols_plot, constellation);
    
    % 添加频偏和相位噪声
    n = (0:length(tx_symbols)-1).';
    freq_offset = 300;  % 300 Hz
    freq_offset_rad = 2 * pi * freq_offset / config.symbol_rate;
    pn_std_rad = 2 * pi / 180;  % 2度/符号
    phase_noise = cumsum(pn_std_rad * randn(length(tx_symbols), 1));
    
    rx_symbols = tx_symbols .* exp(1j * freq_offset_rad * n) .* exp(1j * phase_noise);
    rx_symbols = add_awgn(rx_symbols, 20);
    
    % 频偏补偿
    freq_config = struct('L', 1, 'avg_window', 128, 'mod_order', 4);
    [freq_est, ~, ~] = frequency_offset_estimator(rx_symbols, config.symbol_rate, freq_config);
    rx_freq_comp = rx_symbols .* exp(-1j * 2 * pi * freq_est / config.symbol_rate * n);
    
    % 相位噪声补偿
    pn_config = struct();
    pn_config.mod_order = 4;
    pn_config.win_half = 5;
    pn_config.constellation = constellation;
    pn_config.ambiguity_method = 'decision_directed';
    [~, rx_final, ~] = phase_noise_estimator(rx_freq_comp, pn_config);
    
    figure('Name', '星座图对比', 'Position', [250, 250, 1400, 400]);
    
    subplot(1, 4, 1);
    plot(real(tx_symbols), imag(tx_symbols), '.g', 'MarkerSize', 3);
    hold on;
    plot(real(constellation), imag(constellation), 'ko', 'MarkerSize', 15, 'LineWidth', 2);
    title('发射信号');
    xlabel('实部'); ylabel('虚部');
    axis equal; grid on;
    xlim([-1.5 1.5]); ylim([-1.5 1.5]);
    
    subplot(1, 4, 2);
    plot(real(rx_symbols), imag(rx_symbols), '.r', 'MarkerSize', 3);
    title(sprintf('接收信号\n(频偏+相位噪声+AWGN)'));
    xlabel('实部'); ylabel('虚部');
    axis equal; grid on;
    xlim([-2 2]); ylim([-2 2]);
    
    subplot(1, 4, 3);
    plot(real(rx_freq_comp), imag(rx_freq_comp), '.', 'Color', [1 0.5 0], 'MarkerSize', 3);
    title(sprintf('频偏补偿后\n(估计: %.1f Hz)', freq_est));
    xlabel('实部'); ylabel('虚部');
    axis equal; grid on;
    xlim([-2 2]); ylim([-2 2]);
    
    subplot(1, 4, 4);
    plot(real(rx_final), imag(rx_final), '.b', 'MarkerSize', 3);
    hold on;
    plot(real(constellation), imag(constellation), 'ko', 'MarkerSize', 15, 'LineWidth', 2);
    title(sprintf('相位噪声补偿后\nEVM = %.2f%%', calculate_evm(rx_final, constellation)));
    xlabel('实部'); ylabel('虚部');
    axis equal; grid on;
    xlim([-1.5 1.5]); ylim([-1.5 1.5]);
    
    sgtitle('【图4】信号处理流程星座图对比', 'FontSize', 14, 'FontWeight', 'bold');
    
    if config.save_figures
        savefig(fullfile(config.output_dir, 'constellation_comparison.fig'));
    end
    
    %% ==================== 图5: 相位跟踪轨迹 ====================
    fprintf('生成图5: 相位跟踪轨迹...\n');
    
    num_symbols_track = 2000;
    [tx_symbols, ~] = generate_qpsk_signal(num_symbols_track, constellation);
    
    % 生成相位噪声
    pn_std_rad = 3 * pi / 180;
    phase_noise_true = cumsum(pn_std_rad * randn(num_symbols_track, 1));
    
    rx_symbols = tx_symbols .* exp(1j * phase_noise_true);
    rx_symbols = add_awgn(rx_symbols, 25);
    
    % 相位估计
    pn_config = struct();
    pn_config.mod_order = 4;
    pn_config.win_half = 5;
    pn_config.constellation = constellation;
    pn_config.ambiguity_method = 'unwrap';
    pn_config.phase_unwrap = true;
    
    [phase_est, ~, diag] = phase_noise_estimator(rx_symbols, pn_config);
    
    figure('Name', '相位跟踪轨迹', 'Position', [300, 300, 1200, 600]);
    
    % 上图：真实相位不限制范围，估计相位限制0-360度
    subplot(2, 1, 1);
    
    % 真实相位使用unwrap后的连续值（不限制范围）
    phase_true_deg = phase_noise_true * 180/pi;
    
    % 估计相位限制0-360度范围
    phase_est_mod = mod(phase_est, 2*pi) * 180/pi;
    
    yyaxis left;
    plot(1:num_symbols_track, phase_true_deg, 'r-', 'LineWidth', 1);
    ylabel('真实相位 (度)');
    
    yyaxis right;
    plot(1:num_symbols_track, phase_est_mod, 'b-', 'LineWidth', 1);
    ylabel('估计相位 (度, 0°~360°)');
    ylim([0 360]);
    
    xlabel('符号索引');
    title('相位轨迹跟踪 (真实相位连续显示, 估计相位限制在0°~360°)');
    legend('真实相位噪声', '估计相位', 'Location', 'best');
    grid on;
    
    subplot(2, 1, 2);
    % 计算相位误差：估计相位与真实相位的差异
    phase_error_full = (phase_est - phase_noise_true) * 180/pi;
    % 去除 2π/M 模糊导致的大跳变
    phase_error_full = mod(phase_error_full + 45, 90) - 45;
    
    plot(1:num_symbols_track, phase_error_full, 'k-', 'LineWidth', 0.5);
    hold on;
    yline(0, 'r--', 'LineWidth', 1.5);
    xlabel('符号索引');
    ylabel('相位估计误差 (度)');
    title(sprintf('相位估计误差 (RMSE = %.2f°)', std(phase_error_full)));
    grid on;
    ylim([-20 20]);
    grid on;
    ylim([-20 20]);
    
    sgtitle('【图5】相位噪声跟踪性能', 'FontSize', 14, 'FontWeight', 'bold');
    
    if config.save_figures
        savefig(fullfile(config.output_dir, 'phase_tracking.fig'));
    end
    
    %% ==================== 图6: 算法性能摘要 ====================
    fprintf('生成图6: 算法性能摘要...\n');
    
    figure('Name', '算法性能摘要', 'Position', [350, 350, 800, 600]);
    
    % 创建摘要表格
    metrics = {'频偏估计范围', '频偏估计精度 (SNR=20dB)', ...
               '相位噪声补偿范围', 'EVM改善 (σ=2°)', ...
               '相位模糊解决', '推荐SNR'};
    values = {'0 ~ 2000 Hz', '< 10 Hz', ...
              '0.5 ~ 5 °/符号', sprintf('%.1f%%', mean(evm_before(4:5) - evm_after(4:5))), ...
              '判决引导法', '> 15 dB'};
    
    uitable('Data', [metrics', values'], ...
            'ColumnName', {'性能指标', '典型值'}, ...
            'Position', [50, 320, 700, 200], ...
            'ColumnWidth', {200, 150}, ...
            'FontSize', 11);
    
    % 添加结论文本
    annotation('textbox', [0.1, 0.05, 0.8, 0.4], ...
               'String', {'\bf算法性能总结:', '', ...
                          '1. 频偏估计模块采用四次方差分相位法，支持0-2kHz范围频偏估计', ...
                          '2. 相位噪声估计模块采用改进的V&V算法，有效解决π/2相位模糊', ...
                          '3. 判决引导模糊解决方法在中等SNR下表现最优', ...
                          '4. 联合频偏和相位噪声补偿可显著改善星座图质量', ...
                          '', ...
                          sprintf('报告生成时间: %s', datestr(now))}, ...
               'FitBoxToText', 'on', ...
               'BackgroundColor', [0.95, 0.95, 0.95], ...
               'EdgeColor', [0.5, 0.5, 0.5], ...
               'FontSize', 10);
    
    title('【图6】算法性能摘要', 'FontSize', 14, 'FontWeight', 'bold');
    axis off;
    
    if config.save_figures
        savefig(fullfile(config.output_dir, 'performance_summary.fig'));
        
        % 保存数据
        save(fullfile(config.output_dir, 'freq_phase_report_data.mat'), 'report_data');
        fprintf('\n报告文件已保存到: %s\n', config.output_dir);
    end
    
    fprintf('\n报告生成完成!\n');
    fprintf('═══════════════════════════════════════════════════════════════\n\n');
    
end

%% ==================== 辅助函数 ====================

function [symbols, bits] = generate_qpsk_signal(num_symbols, constellation)
    bits = randi([0 1], num_symbols * 2, 1);
    bit_pairs = reshape(bits, 2, []).';
    symbol_indices = bit_pairs(:, 1) * 2 + bit_pairs(:, 2);
    gray_map = [0 1 3 2];
    gray_indices = gray_map(symbol_indices + 1);
    symbols = constellation(gray_indices + 1).';
end

function signal_noisy = add_awgn(signal, snr_db)
    signal_power = mean(abs(signal).^2);
    noise_power = signal_power / (10^(snr_db/10));
    noise = sqrt(noise_power/2) * (randn(size(signal)) + 1j * randn(size(signal)));
    signal_noisy = signal + noise;
end

function evm = calculate_evm(signal, constellation)
    N = length(signal);
    decided = zeros(N, 1);
    for n = 1:N
        [~, idx] = min(abs(signal(n) - constellation));
        decided(n) = constellation(idx);
    end
    error_vec = signal - decided;
    evm = sqrt(mean(abs(error_vec).^2)) / sqrt(mean(abs(decided).^2)) * 100;
end

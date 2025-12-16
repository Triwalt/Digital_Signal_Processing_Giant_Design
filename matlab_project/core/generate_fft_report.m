% generate_fft_report.m: 生成 my_fft_mix 性能报告图表（用于PPT展示）
%
% 【输出内容】
%   1. 准确性对比图：不同 N 值下的最大绝对误差（对数坐标）
%   2. 性能对比图：my_fft_mix vs MATLAB fft 的执行时间
%   3. 性能比率图：时间比率随 N 的变化
%   4. 综合表格：关键数据汇总
%
% 【使用方法】
%   运行此脚本，将在 reports/ 目录生成：
%   - fft_mix_accuracy.png     准确性对比图
%   - fft_mix_performance.png  性能对比图
%   - fft_mix_summary.png      综合汇总图（单页PPT用）
%   - fft_mix_data.mat         原始数据（供进一步分析）

clear; clc; close all;

fprintf('========================================\n');
fprintf('生成 my_fft_mix 性能报告图表\n');
fprintf('========================================\n\n');

% =====================================================================
% 测试配置
% =====================================================================
% 测试点数：覆盖质数、2的幂、一般合数
test_lengths = [15, 32, 64, 80, 96, 100, 128, 255, 256, 257, 512, 1000, 1024];

% 信号类型
signal_types = {'impulse', 'single_tone', 'white_noise'};
signal_labels = {'冲激', '单频正弦', '白噪声'};

% 性能测试配置
perf_repeats = 20;
perf_warmup = 3;

% 输出目录
output_dir = fullfile(fileparts(mfilename('fullpath')), '..', 'reports');
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

% =====================================================================
% 数据收集
% =====================================================================
fprintf('【阶段1】收集准确性数据...\n');

num_N = length(test_lengths);
num_sig = length(signal_types);

% 存储结果
accuracy_data = struct();
accuracy_data.N = test_lengths;
accuracy_data.signal_types = signal_types;
accuracy_data.max_abs_err = zeros(num_N, num_sig);
accuracy_data.rms_err = zeros(num_N, num_sig);

for i = 1:num_N
    N = test_lengths(i);
    for j = 1:num_sig
        x = generate_test_signal(N, signal_types{j});
        Y_custom = my_fft_mix(x);
        Y_ref = fft(x);
        
        err = abs(Y_custom - Y_ref);
        accuracy_data.max_abs_err(i, j) = max(err);
        accuracy_data.rms_err(i, j) = sqrt(mean(err.^2));
    end
    fprintf('  N=%4d 完成\n', N);
end

fprintf('\n【阶段2】收集性能数据...\n');

perf_data = struct();
perf_data.N = test_lengths;
perf_data.time_custom = zeros(num_N, 1);
perf_data.time_builtin = zeros(num_N, 1);
perf_data.ratio = zeros(num_N, 1);

for i = 1:num_N
    N = test_lengths(i);
    x = randn(1, N) + 1j * randn(1, N);
    
    % Warmup
    for w = 1:perf_warmup
        my_fft_mix(x);
        fft(x);
    end
    
    % 计时 my_fft_mix
    times_custom = zeros(perf_repeats, 1);
    for r = 1:perf_repeats
        tic;
        my_fft_mix(x);
        times_custom(r) = toc;
    end
    
    % 计时 fft
    times_builtin = zeros(perf_repeats, 1);
    for r = 1:perf_repeats
        tic;
        fft(x);
        times_builtin(r) = toc;
    end
    
    perf_data.time_custom(i) = median(times_custom) * 1000;  % ms
    perf_data.time_builtin(i) = median(times_builtin) * 1000;
    perf_data.ratio(i) = perf_data.time_custom(i) / perf_data.time_builtin(i);
    
    fprintf('  N=%4d: my_fft_mix=%.3f ms, fft=%.3f ms, ratio=%.1fx\n', ...
        N, perf_data.time_custom(i), perf_data.time_builtin(i), perf_data.ratio(i));
end

% =====================================================================
% 图1：综合汇总图（单页PPT）
% =====================================================================
fprintf('\n【阶段3】生成图表...\n');

fig_summary = figure('Name', 'FFT Mix 综合报告', 'Color', 'w', ...
    'Position', [100, 100, 1400, 900]);

% 子图布局: 2x2
% (1,1) 准确性对比 - 不同信号类型
% (1,2) 性能对比 - 执行时间
% (2,1) 性能比率
% (2,2) 数据表格

% --- 子图1: 准确性对比 ---
ax1 = subplot(2, 2, 1);
colors = lines(num_sig);
markers = {'o', 's', '^'};

% 将零误差替换为 eps/10 以便在对数坐标显示
plot_err = accuracy_data.max_abs_err;
plot_err(plot_err < eps) = eps / 10;

for j = 1:num_sig
    semilogy(1:num_N, plot_err(:, j), ...
        ['-' markers{j}], 'Color', colors(j,:), 'LineWidth', 1.5, ...
        'MarkerSize', 7, 'MarkerFaceColor', colors(j,:), ...
        'DisplayName', signal_labels{j});
    hold on;
end

% 添加机器精度参考线
yline(eps, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1, 'Label', '机器精度 ε');
hold off;

% 强制设置为对数坐标
set(ax1, 'YScale', 'log');

set(gca, 'XTick', 1:num_N, 'XTickLabel', arrayfun(@num2str, test_lengths, 'UniformOutput', false));
xtickangle(45);
xlabel('FFT 点数 N', 'FontSize', 11);
ylabel('最大绝对误差', 'FontSize', 11);
title('准确性: my\_fft\_mix vs MATLAB fft', 'FontSize', 12, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 9);
grid on;

% 自适应 Y 轴范围：基于数据动态设置
all_err = plot_err(:);
err_min = min(all_err);
err_max = max(all_err);
ylim_low = 10^(floor(log10(err_min)) - 0.5);
ylim_high = 10^(ceil(log10(err_max)) + 0.5);
ylim([ylim_low, ylim_high]);

% --- 子图2: 性能对比 ---
subplot(2, 2, 2);
bar_data = [perf_data.time_custom, perf_data.time_builtin];
b = bar(bar_data, 'grouped');
b(1).FaceColor = [0.2 0.6 0.8];
b(2).FaceColor = [0.9 0.4 0.3];
set(gca, 'XTickLabel', arrayfun(@num2str, test_lengths, 'UniformOutput', false));
xtickangle(45);
xlabel('FFT 点数 N', 'FontSize', 11);
ylabel('执行时间 (ms)', 'FontSize', 11);
title('性能对比: 执行时间', 'FontSize', 12, 'FontWeight', 'bold');
legend({'my\_fft\_mix', 'MATLAB fft'}, 'Location', 'northwest', 'FontSize', 9);
grid on;
set(gca, 'YScale', 'log');

% --- 子图3: 性能比率 ---
subplot(2, 2, 3);

% 区分不同类型的 N
is_power2 = arrayfun(@(n) n == 2^round(log2(n)), test_lengths);
is_prime = arrayfun(@isprime, test_lengths);
is_composite = ~is_power2 & ~is_prime;

hold on;
% 2的幂
idx_p2 = find(is_power2);
bar(idx_p2, perf_data.ratio(is_power2), 0.6, 'FaceColor', [0.3 0.7 0.3], 'EdgeColor', 'none');
% 质数
idx_pr = find(is_prime);
bar(idx_pr, perf_data.ratio(is_prime), 0.6, 'FaceColor', [0.8 0.4 0.4], 'EdgeColor', 'none');
% 一般合数
idx_co = find(is_composite);
bar(idx_co, perf_data.ratio(is_composite), 0.6, 'FaceColor', [0.4 0.5 0.8], 'EdgeColor', 'none');
hold off;

set(gca, 'XTick', 1:num_N, 'XTickLabel', arrayfun(@num2str, test_lengths, 'UniformOutput', false));
xtickangle(45);
xlabel('FFT 点数 N', 'FontSize', 11);
ylabel('时间比率 (my\_fft\_mix / fft)', 'FontSize', 11);
title('性能比率 (按 N 类型着色)', 'FontSize', 12, 'FontWeight', 'bold');
legend({'2的幂', '质数 (Bluestein)', '一般合数'}, 'Location', 'northwest', 'FontSize', 9);
grid on;
yline(1, '--k', 'LineWidth', 1);

% --- 子图4: 数据汇总表格 ---
subplot(2, 2, 4);
axis off;

% 选取代表性的 N 值展示
display_idx = [1, 3, 4, 6, 8, 10, 11, 13];  % 15, 64, 80, 100, 255, 257, 512, 1024
display_idx = display_idx(display_idx <= num_N);

% 准备表格数据
table_data = cell(length(display_idx) + 1, 5);
table_data(1, :) = {'N', '类型', '最大误差', '时间(ms)', '比率'};

for k = 1:length(display_idx)
    i = display_idx(k);
    N = test_lengths(i);
    
    if is_power2(i)
        type_str = '2^n';
    elseif is_prime(i)
        type_str = '质数';
    else
        type_str = '合数';
    end
    
    table_data{k+1, 1} = num2str(N);
    table_data{k+1, 2} = type_str;
    table_data{k+1, 3} = sprintf('%.1e', max(accuracy_data.max_abs_err(i, :)));
    table_data{k+1, 4} = sprintf('%.2f', perf_data.time_custom(i));
    table_data{k+1, 5} = sprintf('%.1fx', perf_data.ratio(i));
end

% 绘制表格
[nrows, ncols] = size(table_data);
col_widths = [0.12, 0.15, 0.25, 0.22, 0.18];
row_height = 0.08;
start_x = 0.08;
start_y = 0.85;

for row = 1:nrows
    for col = 1:ncols
        x = start_x + sum(col_widths(1:col-1));
        y = start_y - (row - 1) * row_height;
        
        if row == 1
            % 表头
            text(x, y, table_data{row, col}, 'FontSize', 10, 'FontWeight', 'bold', ...
                'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle');
        else
            text(x, y, table_data{row, col}, 'FontSize', 9, ...
                'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle');
        end
    end
end

% 表格标题
text(0.5, 0.95, '关键数据汇总', 'FontSize', 12, 'FontWeight', 'bold', ...
    'HorizontalAlignment', 'center');

% 添加说明
text(0.5, 0.08, {'算法特点:', ...
    '• 混合基 Cooley-Tukey: 适用于任意合数 N', ...
    '• Bluestein (Chirp-Z): 用于大质数 N', ...
    '• 误差在机器精度量级 (10^{-12} ~ 10^{-15})'}, ...
    'FontSize', 9, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');

% 整体标题
sgtitle('混合基 FFT (my\_fft\_mix) 性能报告', 'FontSize', 14, 'FontWeight', 'bold');

% 保存图表
saveas(fig_summary, fullfile(output_dir, 'fft_mix_summary.png'));
saveas(fig_summary, fullfile(output_dir, 'fft_mix_summary.fig'));
fprintf('  已保存: fft_mix_summary.png\n');

% =====================================================================
% 图2：单独的准确性图（高清）
% =====================================================================
fig_accuracy = figure('Name', '准确性分析', 'Color', 'w', ...
    'Position', [100, 100, 900, 500]);

% 将零误差替换为 eps/10 以便在对数坐标显示
plot_err_acc = accuracy_data.max_abs_err;
plot_err_acc(plot_err_acc < eps) = eps / 10;

x_base = 1:num_N;
semilogy(x_base, plot_err_acc(:, 1), '-o', 'LineWidth', 2, ...
    'MarkerSize', 8, 'MarkerFaceColor', 'auto', 'DisplayName', signal_labels{1});
hold on;
semilogy(x_base, plot_err_acc(:, 2), '-s', 'LineWidth', 2, ...
    'MarkerSize', 8, 'MarkerFaceColor', 'auto', 'DisplayName', signal_labels{2});
semilogy(x_base, plot_err_acc(:, 3), '-^', 'LineWidth', 2, ...
    'MarkerSize', 8, 'MarkerFaceColor', 'auto', 'DisplayName', signal_labels{3});
yline(eps, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.5, 'Label', '机器精度');
hold off;

% 强制设置为对数坐标
set(gca, 'YScale', 'log');

set(gca, 'XTick', 1:num_N, 'XTickLabel', arrayfun(@num2str, test_lengths, 'UniformOutput', false));
xtickangle(45);
xlabel('FFT 点数 N', 'FontSize', 12);
ylabel('最大绝对误差', 'FontSize', 12);
title('my\_fft\_mix 准确性分析', 'FontSize', 14, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 10);
grid on;

% 自适应 Y 轴范围
all_err_acc = plot_err_acc(:);
err_min = min(all_err_acc);
err_max = max(all_err_acc);
ylim_low = 10^(floor(log10(err_min)) - 0.5);
ylim_high = 10^(ceil(log10(err_max)) + 0.5);
ylim([ylim_low, ylim_high]);

set(gca, 'FontSize', 11);

saveas(fig_accuracy, fullfile(output_dir, 'fft_mix_accuracy.png'));
fprintf('  已保存: fft_mix_accuracy.png\n');

% =====================================================================
% 图3：单独的性能图（高清）
% =====================================================================
fig_perf = figure('Name', '性能分析', 'Color', 'w', ...
    'Position', [100, 100, 900, 500]);

yyaxis left;
loglog(test_lengths, perf_data.time_custom, '-o', 'LineWidth', 2, ...
    'MarkerSize', 8, 'MarkerFaceColor', 'auto');
hold on;
loglog(test_lengths, perf_data.time_builtin, '-s', 'LineWidth', 2, ...
    'MarkerSize', 8, 'MarkerFaceColor', 'auto');
hold off;
ylabel('执行时间 (ms)', 'FontSize', 12);

yyaxis right;
semilogx(test_lengths, perf_data.ratio, '--^', 'LineWidth', 1.5, ...
    'MarkerSize', 6, 'Color', [0.5 0.5 0.5]);
ylabel('时间比率', 'FontSize', 12);

xlabel('FFT 点数 N', 'FontSize', 12);
title('my\_fft\_mix 性能分析', 'FontSize', 14, 'FontWeight', 'bold');
legend({'my\_fft\_mix', 'MATLAB fft', '比率'}, 'Location', 'northwest', 'FontSize', 10);
grid on;
set(gca, 'FontSize', 11);

saveas(fig_perf, fullfile(output_dir, 'fft_mix_performance.png'));
fprintf('  已保存: fft_mix_performance.png\n');

% =====================================================================
% 保存原始数据
% =====================================================================
save(fullfile(output_dir, 'fft_mix_data.mat'), 'accuracy_data', 'perf_data', 'test_lengths');
fprintf('  已保存: fft_mix_data.mat\n');

% =====================================================================
% 输出文本报告
% =====================================================================
fprintf('\n========================================\n');
fprintf('【报告摘要】\n');
fprintf('========================================\n');
fprintf('测试点数范围: %d ~ %d\n', min(test_lengths), max(test_lengths));

% 找到最大误差对应的 N
[max_err, max_idx] = max(accuracy_data.max_abs_err(:));
[row_idx, ~] = ind2sub(size(accuracy_data.max_abs_err), max_idx);
fprintf('最大绝对误差: %.2e (N=%d)\n', max_err, test_lengths(row_idx));
fprintf('最小绝对误差: %.2e\n', min(accuracy_data.max_abs_err(accuracy_data.max_abs_err > 0)));
fprintf('平均性能比率: %.1fx (相对于 MATLAB fft)\n', mean(perf_data.ratio));

[max_ratio, max_ratio_idx] = max(perf_data.ratio);
fprintf('最大性能比率: %.1fx (N=%d)\n', max_ratio, test_lengths(max_ratio_idx));
fprintf('\n输出文件位置: %s\n', output_dir);
fprintf('========================================\n');

% =====================================================================
% 辅助函数
% =====================================================================
function x = generate_test_signal(N, type)
    rng(42);
    n = (0:N-1).';
    
    switch type
        case 'impulse'
            x = zeros(N, 1);
            x(1) = 1;
        case 'single_tone'
            k = max(1, floor(N/8));
            x = exp(1j * 2 * pi * k * n / N);
        case 'white_noise'
            x = randn(N, 1) + 1j * randn(N, 1);
        otherwise
            x = randn(N, 1) + 1j * randn(N, 1);
    end
    x = x(:).';
end

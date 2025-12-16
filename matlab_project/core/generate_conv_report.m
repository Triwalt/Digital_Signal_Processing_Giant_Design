% generate_conv_report.m: 生成快速卷积算法性能报告图表
%
% 【测试对象】
%   1. fast_conv_os.m - 重叠保留法 (Overlap-Save)
%   2. fast_conv_add.m - 重叠相加法 (Overlap-Add)
%   对比基准: MATLAB conv() 函数
%
% 【输出内容】
%   1. 准确性对比图：不同信号长度下的最大绝对误差
%   2. 性能对比图：执行时间 vs 信号长度
%   3. 分块长度影响图：不同 FFT 块长对性能的影响
%   4. 综合报告图（单页PPT用）
%
% 【使用方法】
%   运行此脚本，将在 reports/ 目录生成图表

clear; clc; close all;

fprintf('==========================================\n');
fprintf('生成快速卷积算法性能报告图表\n');
fprintf('==========================================\n\n');

% =====================================================================
% 测试配置
% =====================================================================

% 信号长度测试点
signal_lengths = [1000, 2000, 5000, 10000, 20000, 50000, 100000];

% 滤波器长度测试点
filter_lengths = [32, 64, 128, 256, 512];

% 默认测试参数
default_signal_len = 50000;
default_filter_len = 128;
default_nfft = 1024;  % Overlap-Save 使用
default_L = 512;      % Overlap-Add 使用

% FFT 块长度测试
nfft_values = [256, 512, 1024, 2048, 4096];

% 性能测试配置
perf_repeats = 20;
perf_warmup = 3;

% 输出目录
output_dir = fullfile(fileparts(mfilename('fullpath')), '..', 'reports');
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

% =====================================================================
% 阶段1: 准确性测试
% =====================================================================
fprintf('【阶段1】准确性测试...\n');

accuracy_data = struct();
accuracy_data.signal_lengths = signal_lengths;
accuracy_data.filter_lengths = filter_lengths;
accuracy_data.err_os = zeros(length(signal_lengths), length(filter_lengths));
accuracy_data.err_add = zeros(length(signal_lengths), length(filter_lengths));

rng(42);

for i = 1:length(signal_lengths)
    Nx = signal_lengths(i);
    for j = 1:length(filter_lengths)
        M = filter_lengths(j);
        
        % 生成测试信号
        x = randn(Nx, 1);
        h = randn(M, 1);
        
        % 参考结果
        y_ref = conv(x, h);
        
        % Overlap-Save
        nfft_os = 2^nextpow2(4 * M);  % 典型选择
        try
            y_os = fast_conv_os(x, h, nfft_os);
            accuracy_data.err_os(i, j) = max(abs(y_os(:) - y_ref(:)));
        catch
            accuracy_data.err_os(i, j) = NaN;
        end
        
        % Overlap-Add
        L_add = nfft_os - M + 1;  % 匹配的块长
        try
            y_add = fast_conv_add(x, h, L_add);
            accuracy_data.err_add(i, j) = max(abs(y_add(:) - y_ref(:)));
        catch
            accuracy_data.err_add(i, j) = NaN;
        end
    end
    fprintf('  信号长度 Nx=%6d 完成\n', Nx);
end

% =====================================================================
% 阶段2: 性能测试 - 不同信号长度
% =====================================================================
fprintf('\n【阶段2】性能测试 (变化信号长度)...\n');

perf_data_siglen = struct();
perf_data_siglen.signal_lengths = signal_lengths;
perf_data_siglen.time_os = zeros(length(signal_lengths), 1);
perf_data_siglen.time_add = zeros(length(signal_lengths), 1);
perf_data_siglen.time_conv = zeros(length(signal_lengths), 1);

M = default_filter_len;
h = randn(M, 1);

for i = 1:length(signal_lengths)
    Nx = signal_lengths(i);
    x = randn(Nx, 1);
    
    % Warmup
    for w = 1:perf_warmup
        fast_conv_os(x, h, default_nfft);
        fast_conv_add(x, h, default_L);
        conv(x, h);
    end
    
    % Overlap-Save
    times = zeros(perf_repeats, 1);
    for r = 1:perf_repeats
        tic;
        fast_conv_os(x, h, default_nfft);
        times(r) = toc;
    end
    perf_data_siglen.time_os(i) = median(times) * 1000;
    
    % Overlap-Add
    for r = 1:perf_repeats
        tic;
        fast_conv_add(x, h, default_L);
        times(r) = toc;
    end
    perf_data_siglen.time_add(i) = median(times) * 1000;
    
    % MATLAB conv
    for r = 1:perf_repeats
        tic;
        conv(x, h);
        times(r) = toc;
    end
    perf_data_siglen.time_conv(i) = median(times) * 1000;
    
    fprintf('  Nx=%6d: OS=%.2f ms, Add=%.2f ms, conv=%.2f ms\n', ...
        Nx, perf_data_siglen.time_os(i), perf_data_siglen.time_add(i), ...
        perf_data_siglen.time_conv(i));
end

% =====================================================================
% 阶段3: 性能测试 - 不同滤波器长度
% =====================================================================
fprintf('\n【阶段3】性能测试 (变化滤波器长度)...\n');

perf_data_filtlen = struct();
perf_data_filtlen.filter_lengths = filter_lengths;
perf_data_filtlen.time_os = zeros(length(filter_lengths), 1);
perf_data_filtlen.time_add = zeros(length(filter_lengths), 1);
perf_data_filtlen.time_conv = zeros(length(filter_lengths), 1);

Nx = default_signal_len;
x = randn(Nx, 1);

for j = 1:length(filter_lengths)
    M = filter_lengths(j);
    h = randn(M, 1);
    nfft = 2^nextpow2(4 * M);
    L = nfft - M + 1;
    
    % Warmup
    for w = 1:perf_warmup
        fast_conv_os(x, h, nfft);
        fast_conv_add(x, h, L);
        conv(x, h);
    end
    
    % Overlap-Save
    times = zeros(perf_repeats, 1);
    for r = 1:perf_repeats
        tic;
        fast_conv_os(x, h, nfft);
        times(r) = toc;
    end
    perf_data_filtlen.time_os(j) = median(times) * 1000;
    
    % Overlap-Add
    for r = 1:perf_repeats
        tic;
        fast_conv_add(x, h, L);
        times(r) = toc;
    end
    perf_data_filtlen.time_add(j) = median(times) * 1000;
    
    % MATLAB conv
    for r = 1:perf_repeats
        tic;
        conv(x, h);
        times(r) = toc;
    end
    perf_data_filtlen.time_conv(j) = median(times) * 1000;
    
    fprintf('  M=%4d: OS=%.2f ms, Add=%.2f ms, conv=%.2f ms\n', ...
        M, perf_data_filtlen.time_os(j), perf_data_filtlen.time_add(j), ...
        perf_data_filtlen.time_conv(j));
end

% =====================================================================
% 阶段4: 块长度影响测试
% =====================================================================
fprintf('\n【阶段4】FFT 块长度影响测试...\n');

block_data = struct();
block_data.nfft_values = nfft_values;
block_data.time_os = zeros(length(nfft_values), 1);
block_data.time_add = zeros(length(nfft_values), 1);

Nx = default_signal_len;
M = default_filter_len;
x = randn(Nx, 1);
h = randn(M, 1);

for k = 1:length(nfft_values)
    nfft = nfft_values(k);
    L = nfft - M + 1;
    
    if L < 1
        block_data.time_os(k) = NaN;
        block_data.time_add(k) = NaN;
        continue;
    end
    
    % Warmup
    for w = 1:perf_warmup
        fast_conv_os(x, h, nfft);
        fast_conv_add(x, h, L);
    end
    
    % Overlap-Save
    times = zeros(perf_repeats, 1);
    for r = 1:perf_repeats
        tic;
        fast_conv_os(x, h, nfft);
        times(r) = toc;
    end
    block_data.time_os(k) = median(times) * 1000;
    
    % Overlap-Add
    for r = 1:perf_repeats
        tic;
        fast_conv_add(x, h, L);
        times(r) = toc;
    end
    block_data.time_add(k) = median(times) * 1000;
    
    fprintf('  nfft=%4d (L=%4d): OS=%.2f ms, Add=%.2f ms\n', ...
        nfft, L, block_data.time_os(k), block_data.time_add(k));
end

% =====================================================================
% 图1：综合汇总图（单页PPT）
% =====================================================================
fprintf('\n【阶段5】生成图表...\n');

fig_summary = figure('Name', '快速卷积综合报告', 'Color', 'w', ...
    'Position', [50, 50, 1700, 800]);

% 使用 tiledlayout 并设置紧凑间距
t = tiledlayout(2, 4, 'TileSpacing', 'tight', 'Padding', 'tight');

% --- 子图1: 准确性对比 (选取一个滤波器长度) ---
ax1 = nexttile;
mid_filt_idx = ceil(length(filter_lengths) / 2);

% 处理零误差：将 < eps 的值替换为 eps/10 以便在对数坐标显示
plot_err_os = accuracy_data.err_os(:, mid_filt_idx);
plot_err_add = accuracy_data.err_add(:, mid_filt_idx);
plot_err_os(plot_err_os < eps | isnan(plot_err_os)) = eps / 10;
plot_err_add(plot_err_add < eps | isnan(plot_err_add)) = eps / 10;

x_vals = signal_lengths/1000;
loglog(x_vals, plot_err_os, '-o', ...
    'LineWidth', 2, 'MarkerSize', 7, 'MarkerFaceColor', 'auto', ...
    'DisplayName', 'Overlap-Save');
hold on;
loglog(x_vals, plot_err_add, '-s', ...
    'LineWidth', 2, 'MarkerSize', 7, 'MarkerFaceColor', 'auto', ...
    'DisplayName', 'Overlap-Add');
yline(eps, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1, 'Label', '机器精度');
hold off;

set(ax1, 'XScale', 'log');
set(ax1, 'YScale', 'log');

xlabel('信号长度 (×1000)', 'FontSize', 10);
ylabel('最大绝对误差', 'FontSize', 10);
title(sprintf('准确性 vs MATLAB conv (M=%d)', filter_lengths(mid_filt_idx)), ...
    'FontSize', 11, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 9);
grid on;

% 自适应 Y 轴
all_err = [plot_err_os; plot_err_add];
err_min = min(all_err);
err_max = max(all_err);
ylim([10^(floor(log10(err_min))-0.5), 10^(ceil(log10(err_max))+0.5)]);

% --- 子图2: 性能 vs 信号长度 ---
nexttile;
loglog(signal_lengths/1000, perf_data_siglen.time_os, '-o', 'LineWidth', 2, ...
    'MarkerSize', 7, 'MarkerFaceColor', 'auto', 'DisplayName', 'Overlap-Save');
hold on;
loglog(signal_lengths/1000, perf_data_siglen.time_add, '-s', 'LineWidth', 2, ...
    'MarkerSize', 7, 'MarkerFaceColor', 'auto', 'DisplayName', 'Overlap-Add');
loglog(signal_lengths/1000, perf_data_siglen.time_conv, '-^', 'LineWidth', 2, ...
    'MarkerSize', 7, 'MarkerFaceColor', 'auto', 'DisplayName', 'MATLAB conv');
hold off;

xlabel('信号长度 (×1000)', 'FontSize', 10);
ylabel('执行时间 (ms)', 'FontSize', 10);
title(sprintf('执行时间 vs 信号长度 (M=%d)', default_filter_len), ...
    'FontSize', 11, 'FontWeight', 'bold');
legend('Location', 'northwest', 'FontSize', 9);
grid on;

% --- 子图3: 性能 vs 滤波器长度 (双Y轴) ---
ax3 = nexttile;

% 左Y轴: Overlap-Save 和 Overlap-Add
yyaxis left;
bar_data_left = [perf_data_filtlen.time_os, perf_data_filtlen.time_add];
b_left = bar(bar_data_left, 'grouped');
b_left(1).FaceColor = [0.2 0.6 0.8];
b_left(2).FaceColor = [0.4 0.8 0.4];
ylabel('OS/Add 时间 (ms)', 'FontSize', 9);
set(gca, 'YColor', 'k');

% 右Y轴: MATLAB conv
yyaxis right;
plot(1:length(filter_lengths), perf_data_filtlen.time_conv, '-o', ...
    'LineWidth', 2, 'MarkerSize', 6, 'Color', [0.9 0.4 0.3], ...
    'MarkerFaceColor', [0.9 0.4 0.3]);
ylabel('MATLAB conv 时间 (ms)', 'FontSize', 9, 'Color', [0.9 0.4 0.3]);
set(gca, 'YColor', [0.9 0.4 0.3]);

set(gca, 'XTick', 1:length(filter_lengths));
set(gca, 'XTickLabel', arrayfun(@num2str, filter_lengths, 'UniformOutput', false));
xlabel('滤波器长度 M', 'FontSize', 10);
title(sprintf('执行时间 vs 滤波器长度 (Nx=%dk)', default_signal_len/1000), ...
    'FontSize', 11, 'FontWeight', 'bold');
legend({'Overlap-Save', 'Overlap-Add', 'MATLAB conv'}, 'Location', 'northeast', 'FontSize', 8);
grid on;

% --- 子图4: 块长度影响 ---
nexttile;
plot(nfft_values, block_data.time_os, '-o', 'LineWidth', 2, ...
    'MarkerSize', 8, 'MarkerFaceColor', 'auto', 'DisplayName', 'Overlap-Save');
hold on;
plot(nfft_values, block_data.time_add, '-s', 'LineWidth', 2, ...
    'MarkerSize', 8, 'MarkerFaceColor', 'auto', 'DisplayName', 'Overlap-Add');
hold off;

set(gca, 'XScale', 'log');
set(gca, 'XTick', nfft_values);
set(gca, 'XTickLabel', arrayfun(@num2str, nfft_values, 'UniformOutput', false));
xlabel('FFT 块长度 (nfft)', 'FontSize', 10);
ylabel('执行时间 (ms)', 'FontSize', 10);
title(sprintf('块长度影响 (Nx=%dk, M=%d)', default_signal_len/1000, default_filter_len), ...
    'FontSize', 11, 'FontWeight', 'bold');
legend('Location', 'northeast', 'FontSize', 9);
grid on;

% --- 子图5: 加速比 (vs MATLAB conv) ---
nexttile;
speedup_os = perf_data_siglen.time_conv ./ perf_data_siglen.time_os;
speedup_add = perf_data_siglen.time_conv ./ perf_data_siglen.time_add;

semilogx(signal_lengths/1000, speedup_os, '-o', 'LineWidth', 2, ...
    'MarkerSize', 8, 'MarkerFaceColor', 'auto', 'DisplayName', 'Overlap-Save');
hold on;
semilogx(signal_lengths/1000, speedup_add, '-s', 'LineWidth', 2, ...
    'MarkerSize', 8, 'MarkerFaceColor', 'auto', 'DisplayName', 'Overlap-Add');
yline(1, '--k', 'LineWidth', 1, 'Label', '相等');
hold off;

xlabel('信号长度 (×1000)', 'FontSize', 10);
ylabel('加速比 (conv/算法)', 'FontSize', 10);
title('相对 MATLAB conv 的加速比', 'FontSize', 11, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 9);
grid on;

% 自适应 Y 轴范围
all_speedup = [speedup_os; speedup_add];
speedup_min = min(all_speedup);
speedup_max = max(all_speedup);
y_margin = (speedup_max - speedup_min) * 0.1;
if y_margin < 0.01
    y_margin = speedup_max * 0.1;
end
ylim([max(0, speedup_min - y_margin), speedup_max + y_margin]);

% --- 子图6: Overlap-Save 热图 ---
ax6a = nexttile;
log_err_os = log10(accuracy_data.err_os + eps);
log_err_add = log10(accuracy_data.err_add + eps);
cmin = min([log_err_os(:); log_err_add(:)]);
cmax = max([log_err_os(:); log_err_add(:)]);

imagesc(log_err_os);
caxis([cmin, cmax]);
colormap(ax6a, flipud(hot));
set(gca, 'XTick', 1:length(filter_lengths));
set(gca, 'XTickLabel', arrayfun(@num2str, filter_lengths, 'UniformOutput', false));
set(gca, 'YTick', 1:length(signal_lengths));
set(gca, 'YTickLabel', arrayfun(@(x) sprintf('%dk', x/1000), signal_lengths, 'UniformOutput', false));
xlabel('滤波器长度 M', 'FontSize', 10);
ylabel('信号长度', 'FontSize', 10);
title('误差热图: Overlap-Save', 'FontSize', 11, 'FontWeight', 'bold');

% --- 子图7: Overlap-Add 热图 ---
ax6b = nexttile;
imagesc(log_err_add);
caxis([cmin, cmax]);
colormap(ax6b, flipud(hot));
set(gca, 'XTick', 1:length(filter_lengths));
set(gca, 'XTickLabel', arrayfun(@num2str, filter_lengths, 'UniformOutput', false));
set(gca, 'YTick', 1:length(signal_lengths));
set(gca, 'YTickLabel', arrayfun(@(x) sprintf('%dk', x/1000), signal_lengths, 'UniformOutput', false));
xlabel('滤波器长度 M', 'FontSize', 10);
ylabel('信号长度', 'FontSize', 10);
title('误差热图: Overlap-Add', 'FontSize', 11, 'FontWeight', 'bold');

% --- 子图8: 数据汇总表格 ---
ax_table = nexttile;
axis off;

% 添加共享颜色条
cb = colorbar(ax6b, 'eastoutside');
cb.Label.String = 'log_{10}(误差)';
cb.FontSize = 9;

% 汇总信息文本
text_content = {
    sprintf('【测试配置】'), ...
    sprintf('信号长度: %d ~ %d', min(signal_lengths), max(signal_lengths)), ...
    sprintf('滤波器长度: %d ~ %d', min(filter_lengths), max(filter_lengths)), ...
    '', ...
    sprintf('【准确性 (最大误差)】'), ...
    sprintf('OS: %.1e ~ %.1e', min(accuracy_data.err_os(:)), max(accuracy_data.err_os(:))), ...
    sprintf('Add: %.1e ~ %.1e', min(accuracy_data.err_add(:)), max(accuracy_data.err_add(:))), ...
    '', ...
    sprintf('【平均加速比 vs conv】'), ...
    sprintf('OS: %.3fx  Add: %.3fx', mean(speedup_os), mean(speedup_add)), ...
    '', ...
    sprintf('【M↑ 时间↓ 原因】'), ...
    sprintf('nfft=4M, L=nfft-M+1'), ...
    sprintf('M大→块数少→循环开销低')
};

text(0.05, 0.95, text_content, 'FontSize', 8, 'VerticalAlignment', 'top', ...
    'FontName', 'FixedWidth', 'Interpreter', 'none');
title('汇总信息', 'FontSize', 11, 'FontWeight', 'bold');

% 整体标题
title(t, '快速卷积算法性能报告: Overlap-Save vs Overlap-Add', ...
    'FontSize', 14, 'FontWeight', 'bold');

% 保存
saveas(fig_summary, fullfile(output_dir, 'conv_summary.png'));
saveas(fig_summary, fullfile(output_dir, 'conv_summary.fig'));
fprintf('  已保存: conv_summary.png\n');

% =====================================================================
% 图2：单独的准确性详细图
% =====================================================================
fig_acc = figure('Name', '快速卷积准确性', 'Color', 'w', ...
    'Position', [100, 100, 1200, 500]);

% Overlap-Save 热图
subplot(1, 2, 1);
log_err_os = log10(accuracy_data.err_os + eps);
imagesc(log_err_os);
cb = colorbar;
cb.Label.String = 'log_{10}(误差)';
colormap(flipud(hot));
set(gca, 'XTick', 1:length(filter_lengths));
set(gca, 'XTickLabel', arrayfun(@num2str, filter_lengths, 'UniformOutput', false));
set(gca, 'YTick', 1:length(signal_lengths));
set(gca, 'YTickLabel', arrayfun(@(x) sprintf('%dk', x/1000), signal_lengths, 'UniformOutput', false));
xlabel('滤波器长度 M', 'FontSize', 11);
ylabel('信号长度', 'FontSize', 11);
title('Overlap-Save 准确性', 'FontSize', 12, 'FontWeight', 'bold');

% Overlap-Add 热图
subplot(1, 2, 2);
log_err_add = log10(accuracy_data.err_add + eps);
imagesc(log_err_add);
cb = colorbar;
cb.Label.String = 'log_{10}(误差)';
colormap(flipud(hot));
set(gca, 'XTick', 1:length(filter_lengths));
set(gca, 'XTickLabel', arrayfun(@num2str, filter_lengths, 'UniformOutput', false));
set(gca, 'YTick', 1:length(signal_lengths));
set(gca, 'YTickLabel', arrayfun(@(x) sprintf('%dk', x/1000), signal_lengths, 'UniformOutput', false));
xlabel('滤波器长度 M', 'FontSize', 11);
ylabel('信号长度', 'FontSize', 11);
title('Overlap-Add 准确性', 'FontSize', 12, 'FontWeight', 'bold');

sgtitle('快速卷积算法准确性分析', 'FontSize', 14, 'FontWeight', 'bold');

saveas(fig_acc, fullfile(output_dir, 'conv_accuracy.png'));
fprintf('  已保存: conv_accuracy.png\n');

% =====================================================================
% 保存数据
% =====================================================================
save(fullfile(output_dir, 'conv_data.mat'), ...
    'accuracy_data', 'perf_data_siglen', 'perf_data_filtlen', 'block_data', ...
    'signal_lengths', 'filter_lengths', 'nfft_values', ...
    'default_signal_len', 'default_filter_len');
fprintf('  已保存: conv_data.mat\n');

% =====================================================================
% 报告摘要
% =====================================================================
fprintf('\n==========================================\n');
fprintf('【报告摘要】\n');
fprintf('==========================================\n');
fprintf('测试信号长度范围: %d ~ %d\n', min(signal_lengths), max(signal_lengths));
fprintf('测试滤波器长度范围: %d ~ %d\n', min(filter_lengths), max(filter_lengths));

fprintf('\n准确性 (最大绝对误差):\n');
fprintf('  Overlap-Save: %.2e ~ %.2e\n', min(accuracy_data.err_os(:)), max(accuracy_data.err_os(:)));
fprintf('  Overlap-Add:  %.2e ~ %.2e\n', min(accuracy_data.err_add(:)), max(accuracy_data.err_add(:)));

fprintf('\n性能 (Nx=%dk, M=%d):\n', default_signal_len/1000, default_filter_len);
mid_idx = ceil(length(signal_lengths) / 2);
fprintf('  Overlap-Save: %.2f ms\n', perf_data_siglen.time_os(mid_idx));
fprintf('  Overlap-Add:  %.2f ms\n', perf_data_siglen.time_add(mid_idx));
fprintf('  MATLAB conv:  %.2f ms\n', perf_data_siglen.time_conv(mid_idx));

fprintf('\n平均加速比 (vs MATLAB conv):\n');
fprintf('  Overlap-Save: %.2fx\n', mean(speedup_os));
fprintf('  Overlap-Add:  %.2fx\n', mean(speedup_add));

fprintf('\n输出文件位置: %s\n', output_dir);
fprintf('==========================================\n');

%% generate_cma_comparison_report.m: CMA时域/频域实现对比报告
%
% 功能说明:
%   对比 CMA_homework.m (频域FFT实现) 与 CMA_reference.m (时域实现)
%   在处理 TD_TRdata.mat 数据时的:
%   1. 计算精度 (输出误差)
%   2. 运行速度 (执行时间)
%   3. 收敛特性 (误差曲线)
%   4. 星座图质量
%
% 输出:
%   - 综合报告图 (PNG)
%   - 数据文件 (MAT)
%
% 作者: DSP课程设计项目
% 日期: 2025年12月

function report_data = generate_cma_comparison_report(config)

    clc;
    close all;
    
    %% 路径设置
    scriptDir = fileparts(mfilename('fullpath'));
    if isempty(scriptDir)
        scriptDir = pwd;
    end
    
    % 添加必要路径
    coreDir = fullfile(scriptDir);
    hwDir = fullfile(scriptDir, '..', 'homework1');
    addpath(coreDir);
    addpath(hwDir);
    
    fprintf('╔══════════════════════════════════════════════════════════════╗\n');
    fprintf('║     CMA 时域/频域实现 精度与速度对比报告生成器              ║\n');
    fprintf('╚══════════════════════════════════════════════════════════════╝\n\n');
    
    %% 默认配置
    if nargin < 1 || isempty(config)
        config = struct();
    end
    
    if ~isfield(config, 'save_figures')
        config.save_figures = true;
    end
    if ~isfield(config, 'output_dir')
        config.output_dir = fullfile(scriptDir, '..', 'reports');
    end
    if ~isfield(config, 'num_timing_runs')
        config.num_timing_runs = 3;  % 时间测量重复次数
    end
    
    % 确保输出目录存在
    if config.save_figures && ~exist(config.output_dir, 'dir')
        mkdir(config.output_dir);
    end
    
    % 报告数据
    report_data = struct();
    report_data.config = config;
    report_data.timestamp = datetime('now');
    
    %% ==================== 加载数据 ====================
    fprintf('加载TD_TRdata.mat数据...\n');
    
    dataFile = fullfile(hwDir, 'TD_TRdata.mat');
    if ~exist(dataFile, 'file')
        error('找不到数据文件: %s', dataFile);
    end
    
    % 保存当前目录
    originalDir = pwd;
    
    %% ==================== 运行时域CMA (Reference) ====================
    fprintf('\n[1/4] 运行时域CMA (CMA_reference.m)...\n');
    
    cd(hwDir);
    
    % 预设配置变量，防止 CMA_homework.m 的默认配置影响后续运行
    enable_vv_phase_lock = false;
    enable_plot_output = false;
    
    % 时间测量
    time_ref = zeros(config.num_timing_runs, 1);
    for run_idx = 1:config.num_timing_runs
        % 清除工作区变量（保留必要变量和配置开关）
        clearvars -except config report_data hwDir coreDir originalDir scriptDir time_ref run_idx dataFile ...
                  enable_vv_phase_lock enable_plot_output
        
        tic;
        
        % 运行时域CMA
        run('CMA_reference.m');
        
        time_ref(run_idx) = toc;
        fprintf('  运行 %d/%d: %.3f 秒\n', run_idx, config.num_timing_runs, time_ref(run_idx));
    end
    
    % 保存时域结果
    ref_SigX_out = SigX_out;
    ref_SigY_out = SigY_out;
    % 注意: CMA_reference.m 会覆盖 CMAdataX/Y 为尾部数据，因此从 SigX_out 重新生成完整数据
    ref_CMAdataX = SigX_out.';
    ref_CMAdataY = SigY_out.';
    ref_errx_vec = errx_vec;
    ref_erry_vec = erry_vec;
    ref_xx = xx;
    ref_yy = yy;
    ref_xy = xy;
    ref_yx = yx;
    
    close all;  % 关闭CMA_reference生成的图
    
    %% ==================== 运行频域CMA (Homework) ====================
    fprintf('\n[2/4] 运行频域CMA (CMA_homework.m)...\n');
    
    % 时间测量
    time_hw = zeros(config.num_timing_runs, 1);
    for run_idx = 1:config.num_timing_runs
        % 清除工作区变量（保留必要变量）
        clearvars -except config report_data hwDir coreDir originalDir scriptDir ...
                  time_ref time_hw run_idx dataFile ...
                  ref_SigX_out ref_SigY_out ref_CMAdataX ref_CMAdataY ...
                  ref_errx_vec ref_erry_vec ref_xx ref_yy ref_xy ref_yx ...
                  enable_vv_phase_lock enable_plot_output
        
        % 配置开关 (必须在clearvars之后、run之前设置，确保CMA_homework.m使用正确的配置)
        enable_vv_phase_lock = false;  % 关闭相位校正以公平比较
        enable_plot_output = false;
        
        tic;
        
        % 运行频域CMA
        run('CMA_homework.m');
        
        time_hw(run_idx) = toc;
        fprintf('  运行 %d/%d: %.3f 秒\n', run_idx, config.num_timing_runs, time_hw(run_idx));
    end
    
    % 保存频域结果
    hw_SigX_out = SigX_out;
    hw_SigY_out = SigY_out;
    % 保持与时域结果一致的数据来源
    hw_CMAdataX = SigX_out.';
    hw_CMAdataY = SigY_out.';
    hw_errx_vec = errx_vec;
    hw_erry_vec = erry_vec;
    hw_xx = xx;
    hw_yy = yy;
    hw_xy = xy;
    hw_yx = yx;
    
    close all;  % 关闭CMA_homework生成的图
    
    % 返回原目录
    cd(originalDir);
    
    %% ==================== 精度分析 ====================
    fprintf('\n[3/4] 分析精度差异...\n');
    
    % 确保长度一致
    len = min(length(ref_SigX_out), length(hw_SigX_out));
    ref_X = ref_SigX_out(1:len);
    ref_Y = ref_SigY_out(1:len);
    hw_X = hw_SigX_out(1:len);
    hw_Y = hw_SigY_out(1:len);
    
    % 计算误差
    err_X = abs(ref_X - hw_X);
    err_Y = abs(ref_Y - hw_Y);
    
    % 相对误差 (避免除零)
    ref_mag_X = abs(ref_X);
    ref_mag_Y = abs(ref_Y);
    rel_err_X = err_X ./ max(ref_mag_X, 1e-10);
    rel_err_Y = err_Y ./ max(ref_mag_Y, 1e-10);
    
    % 误差统计
    accuracy = struct();
    accuracy.max_abs_err_X = max(err_X);
    accuracy.max_abs_err_Y = max(err_Y);
    accuracy.mean_abs_err_X = mean(err_X);
    accuracy.mean_abs_err_Y = mean(err_Y);
    accuracy.mse_X = mean(err_X.^2);
    accuracy.mse_Y = mean(err_Y.^2);
    accuracy.max_rel_err_X = max(rel_err_X);
    accuracy.max_rel_err_Y = max(rel_err_Y);
    accuracy.mean_rel_err_X = mean(rel_err_X);
    accuracy.mean_rel_err_Y = mean(rel_err_Y);
    
    % 尾部误差 (稳态)
    tail_len = min(10000, len);
    tail_err_X = err_X(end-tail_len+1:end);
    tail_err_Y = err_Y(end-tail_len+1:end);
    accuracy.tail_max_err_X = max(tail_err_X);
    accuracy.tail_max_err_Y = max(tail_err_Y);
    accuracy.tail_mse_X = mean(tail_err_X.^2);
    accuracy.tail_mse_Y = mean(tail_err_Y.^2);
    
    % 抽头差异
    tap_diff_xx = max(abs(ref_xx(:,end) - hw_xx(:,end)));
    tap_diff_yy = max(abs(ref_yy(:,end) - hw_yy(:,end)));
    tap_diff_xy = max(abs(ref_xy(:,end) - hw_xy(:,end)));
    tap_diff_yx = max(abs(ref_yx(:,end) - hw_yx(:,end)));
    accuracy.max_tap_diff = max([tap_diff_xx, tap_diff_yy, tap_diff_xy, tap_diff_yx]);
    
    fprintf('  最大绝对误差 X: %.2e\n', accuracy.max_abs_err_X);
    fprintf('  最大绝对误差 Y: %.2e\n', accuracy.max_abs_err_Y);
    fprintf('  平均绝对误差 X: %.2e\n', accuracy.mean_abs_err_X);
    fprintf('  平均绝对误差 Y: %.2e\n', accuracy.mean_abs_err_Y);
    fprintf('  最大抽头差异: %.2e\n', accuracy.max_tap_diff);
    
    %% ==================== 速度分析 ====================
    fprintf('\n[4/4] 分析速度差异...\n');
    
    timing = struct();
    timing.ref_times = time_ref;
    timing.hw_times = time_hw;
    timing.ref_mean = mean(time_ref);
    timing.ref_std = std(time_ref);
    timing.hw_mean = mean(time_hw);
    timing.hw_std = std(time_hw);
    timing.speedup = timing.ref_mean / timing.hw_mean;
    
    fprintf('  时域CMA平均时间: %.3f ± %.3f 秒\n', timing.ref_mean, timing.ref_std);
    fprintf('  频域CMA平均时间: %.3f ± %.3f 秒\n', timing.hw_mean, timing.hw_std);
    fprintf('  加速比: %.2fx\n', timing.speedup);
    
    %% ==================== 生成综合报告图 ====================
    fprintf('\n生成综合报告图...\n');
    
    figure('Name', 'CMA时域/频域对比报告', 'Position', [50, 50, 1600, 900]);
    
    % 第一行: 星座图对比
    subplot(3, 4, 1);
    plot(ref_CMAdataX, '.r', 'MarkerSize', 1);
    title('时域CMA - X极化');
    xlabel('实部'); ylabel('虚部');
    axis equal; grid on;
    xlim([-2.5 2.5]); ylim([-2.5 2.5]);
    
    subplot(3, 4, 2);
    plot(hw_CMAdataX, '.b', 'MarkerSize', 1);
    title('频域CMA - X极化');
    xlabel('实部'); ylabel('虚部');
    axis equal; grid on;
    xlim([-2.5 2.5]); ylim([-2.5 2.5]);
    
    subplot(3, 4, 3);
    plot(ref_CMAdataY, '.r', 'MarkerSize', 1);
    title('时域CMA - Y极化');
    xlabel('实部'); ylabel('虚部');
    axis equal; grid on;
    xlim([-2.5 2.5]); ylim([-2.5 2.5]);
    
    subplot(3, 4, 4);
    plot(hw_CMAdataY, '.b', 'MarkerSize', 1);
    title('频域CMA - Y极化');
    xlabel('实部'); ylabel('虚部');
    axis equal; grid on;
    xlim([-2.5 2.5]); ylim([-2.5 2.5]);
    
    % 第二行: 误差分析
    subplot(3, 4, 5);
    semilogy(err_X, 'r-', 'LineWidth', 0.5);
    xlabel('样本索引'); ylabel('绝对误差 (log)');
    title('X极化输出误差');
    grid on;
    
    subplot(3, 4, 6);
    semilogy(err_Y, 'b-', 'LineWidth', 0.5);
    xlabel('样本索引'); ylabel('绝对误差 (log)');
    title('Y极化输出误差');
    grid on;
    
    subplot(3, 4, 7);
    % 尾部100点对比
    tail_100 = 100;
    idx_tail = (len-tail_100+1):len;
    plot(1:tail_100, real(ref_X(idx_tail)), 'r-', 'LineWidth', 1.5);
    hold on;
    plot(1:tail_100, real(hw_X(idx_tail)), 'b--', 'LineWidth', 1.5);
    xlabel('样本索引 (后100点)'); ylabel('实部');
    title('X极化尾部对比 (实部)');
    legend('时域', '频域', 'Location', 'best');
    grid on;
    
    subplot(3, 4, 8);
    plot(1:tail_100, imag(ref_X(idx_tail)), 'r-', 'LineWidth', 1.5);
    hold on;
    plot(1:tail_100, imag(hw_X(idx_tail)), 'b--', 'LineWidth', 1.5);
    xlabel('样本索引 (后100点)'); ylabel('虚部');
    title('X极化尾部对比 (虚部)');
    legend('时域', '频域', 'Location', 'best');
    grid on;
    
    % 第三行: 收敛曲线和性能统计
    subplot(3, 4, 9);
    ref_err_smooth = movmean(abs(ref_errx_vec), 1000);
    hw_err_smooth = movmean(abs(hw_errx_vec), 1000);
    plot(ref_err_smooth, 'r-', 'LineWidth', 1);
    hold on;
    plot(hw_err_smooth, 'b--', 'LineWidth', 1);
    xlabel('样本索引'); ylabel('|误差|');
    title('CMA收敛曲线 (X极化)');
    legend('时域', '频域', 'Location', 'northeast');
    grid on;
    
    subplot(3, 4, 10);
    ref_erry_smooth = movmean(abs(ref_erry_vec), 1000);
    hw_erry_smooth = movmean(abs(hw_erry_vec), 1000);
    plot(ref_erry_smooth, 'r-', 'LineWidth', 1);
    hold on;
    plot(hw_erry_smooth, 'b--', 'LineWidth', 1);
    xlabel('样本索引'); ylabel('|误差|');
    title('CMA收敛曲线 (Y极化)');
    legend('时域', '频域', 'Location', 'northeast');
    grid on;
    
    % 速度对比柱状图
    subplot(3, 4, 11);
    bar_data = [timing.ref_mean, timing.hw_mean];
    bar_err = [timing.ref_std, timing.hw_std];
    b = bar(bar_data);
    hold on;
    errorbar(1:2, bar_data, bar_err, 'k.', 'LineWidth', 1.5);
    set(gca, 'XTickLabel', {'时域CMA', '频域CMA'});
    ylabel('运行时间 (秒)');
    title(sprintf('速度对比 (加速比: %.2fx)', timing.speedup));
    grid on;
    
    % 添加数值标签
    for i = 1:2
        text(i, bar_data(i) + bar_err(i) + 0.1, sprintf('%.2fs', bar_data(i)), ...
             'HorizontalAlignment', 'center', 'FontSize', 10);
    end
    
    % 性能摘要
    subplot(3, 4, 12);
    axis off;
    text(0.05, 0.95, '【精度与性能摘要】', 'FontSize', 12, 'FontWeight', 'bold');
    text(0.05, 0.82, sprintf('数据长度: %d 样本', len), 'FontSize', 10);
    text(0.05, 0.70, sprintf('最大绝对误差: %.2e', max(accuracy.max_abs_err_X, accuracy.max_abs_err_Y)), 'FontSize', 10);
    text(0.05, 0.58, sprintf('平均绝对误差: %.2e', mean([accuracy.mean_abs_err_X, accuracy.mean_abs_err_Y])), 'FontSize', 10);
    text(0.05, 0.46, sprintf('尾部MSE: %.2e', mean([accuracy.tail_mse_X, accuracy.tail_mse_Y])), 'FontSize', 10);
    text(0.05, 0.34, sprintf('最大抽头差异: %.2e', accuracy.max_tap_diff), 'FontSize', 10);
    text(0.05, 0.22, sprintf('时域运行时间: %.3f s', timing.ref_mean), 'FontSize', 10);
    text(0.05, 0.10, sprintf('频域运行时间: %.3f s', timing.hw_mean), 'FontSize', 10);
    text(0.05, -0.02, sprintf('加速比: %.2fx', timing.speedup), 'FontSize', 10, 'FontWeight', 'bold', 'Color', 'b');
    
    % 判断是否通过
    if accuracy.max_abs_err_X < 1e-10 && accuracy.max_abs_err_Y < 1e-10
        result_str = '✓ 精度验证通过';
        result_color = [0 0.6 0];
    else
        result_str = '△ 存在微小差异';
        result_color = [0.8 0.5 0];
    end
    text(0.05, -0.15, result_str, 'FontSize', 11, 'FontWeight', 'bold', 'Color', result_color);
    
    sgtitle('【综合报告】CMA均衡器 时域/频域实现对比', 'FontSize', 14, 'FontWeight', 'bold');
    
    if config.save_figures
        savefig(fullfile(config.output_dir, 'cma_comparison_report.fig'));
        print(fullfile(config.output_dir, 'cma_comparison_report.png'), '-dpng', '-r150');
        fprintf('  保存到: %s\n', fullfile(config.output_dir, 'cma_comparison_report.png'));
    end
    
    %% ==================== 生成详细误差分析图 ====================
    fprintf('生成详细误差分析图...\n');
    
    figure('Name', '误差详细分析', 'Position', [100, 100, 1200, 800]);
    
    % 误差分布直方图
    subplot(2, 3, 1);
    histogram(log10(err_X + 1e-20), 50, 'FaceColor', 'r', 'FaceAlpha', 0.7);
    xlabel('log_{10}(误差)'); ylabel('计数');
    title('X极化误差分布 (对数)');
    grid on;
    
    subplot(2, 3, 2);
    histogram(log10(err_Y + 1e-20), 50, 'FaceColor', 'b', 'FaceAlpha', 0.7);
    xlabel('log_{10}(误差)'); ylabel('计数');
    title('Y极化误差分布 (对数)');
    grid on;
    
    % 误差随时间变化
    subplot(2, 3, 3);
    window_size = 1000;
    err_X_moving = movmean(err_X, window_size);
    err_Y_moving = movmean(err_Y, window_size);
    semilogy(err_X_moving, 'r-', 'LineWidth', 1);
    hold on;
    semilogy(err_Y_moving, 'b-', 'LineWidth', 1);
    xlabel('样本索引'); ylabel('滑动平均误差 (log)');
    title(sprintf('误差趋势 (窗口=%d)', window_size));
    legend('X极化', 'Y极化');
    grid on;
    
    % 抽头对比
    subplot(2, 3, 4);
    tap_idx = 1:length(ref_xx(:,end));
    plot(tap_idx, abs(ref_xx(:,end)), 'r-o', 'LineWidth', 1, 'MarkerSize', 4);
    hold on;
    plot(tap_idx, abs(hw_xx(:,end)), 'b--s', 'LineWidth', 1, 'MarkerSize', 4);
    xlabel('抽头索引'); ylabel('|抽头值|');
    title('XX抽头对比');
    legend('时域', '频域', 'Location', 'best');
    grid on;
    
    subplot(2, 3, 5);
    plot(tap_idx, abs(ref_yy(:,end)), 'r-o', 'LineWidth', 1, 'MarkerSize', 4);
    hold on;
    plot(tap_idx, abs(hw_yy(:,end)), 'b--s', 'LineWidth', 1, 'MarkerSize', 4);
    xlabel('抽头索引'); ylabel('|抽头值|');
    title('YY抽头对比');
    legend('时域', '频域', 'Location', 'best');
    grid on;
    
    % 运行时间详情
    subplot(2, 3, 6);
    x_runs = 1:config.num_timing_runs;
    plot(x_runs, time_ref, 'ro-', 'LineWidth', 2, 'MarkerSize', 8);
    hold on;
    plot(x_runs, time_hw, 'bs-', 'LineWidth', 2, 'MarkerSize', 8);
    xlabel('运行次数'); ylabel('时间 (秒)');
    title('多次运行时间');
    legend('时域CMA', '频域CMA', 'Location', 'best');
    grid on;
    xlim([0.5, config.num_timing_runs + 0.5]);
    
    sgtitle('【详细分析】误差分布与抽头对比', 'FontSize', 14, 'FontWeight', 'bold');
    
    if config.save_figures
        savefig(fullfile(config.output_dir, 'cma_comparison_detail.fig'));
        print(fullfile(config.output_dir, 'cma_comparison_detail.png'), '-dpng', '-r150');
        fprintf('  保存到: %s\n', fullfile(config.output_dir, 'cma_comparison_detail.png'));
    end
    
    %% ==================== 保存报告数据 ====================
    report_data.accuracy = accuracy;
    report_data.timing = timing;
    report_data.data_length = len;
    report_data.err_X = err_X;
    report_data.err_Y = err_Y;
    report_data.ref_taps = struct('xx', ref_xx(:,end), 'yy', ref_yy(:,end), ...
                                   'xy', ref_xy(:,end), 'yx', ref_yx(:,end));
    report_data.hw_taps = struct('xx', hw_xx(:,end), 'yy', hw_yy(:,end), ...
                                  'xy', hw_xy(:,end), 'yx', hw_yx(:,end));
    
    if config.save_figures
        save(fullfile(config.output_dir, 'cma_comparison_data.mat'), 'report_data');
        fprintf('  数据保存到: %s\n', fullfile(config.output_dir, 'cma_comparison_data.mat'));
    end
    
    %% ==================== 打印最终报告 ====================
    fprintf('\n');
    fprintf('═══════════════════════════════════════════════════════════════\n');
    fprintf('                    CMA对比测试报告                            \n');
    fprintf('═══════════════════════════════════════════════════════════════\n');
    fprintf('数据长度: %d 样本\n', len);
    fprintf('───────────────────────────────────────────────────────────────\n');
    fprintf('【精度分析】\n');
    fprintf('  X极化 - 最大误差: %.2e, 平均误差: %.2e\n', accuracy.max_abs_err_X, accuracy.mean_abs_err_X);
    fprintf('  Y极化 - 最大误差: %.2e, 平均误差: %.2e\n', accuracy.max_abs_err_Y, accuracy.mean_abs_err_Y);
    fprintf('  尾部MSE: X=%.2e, Y=%.2e\n', accuracy.tail_mse_X, accuracy.tail_mse_Y);
    fprintf('  最大抽头差异: %.2e\n', accuracy.max_tap_diff);
    fprintf('───────────────────────────────────────────────────────────────\n');
    fprintf('【速度分析】\n');
    fprintf('  时域CMA: %.3f ± %.3f 秒\n', timing.ref_mean, timing.ref_std);
    fprintf('  频域CMA: %.3f ± %.3f 秒\n', timing.hw_mean, timing.hw_std);
    fprintf('  加速比: %.2fx\n', timing.speedup);
    fprintf('───────────────────────────────────────────────────────────────\n');
    
    if accuracy.max_abs_err_X < 1e-10 && accuracy.max_abs_err_Y < 1e-10
        fprintf('【结论】频域CMA与时域CMA输出完全一致，验证通过！\n');
    elseif accuracy.max_abs_err_X < 1e-6 && accuracy.max_abs_err_Y < 1e-6
        fprintf('【结论】频域CMA与时域CMA存在微小数值差异（<1e-6），精度可接受\n');
    else
        fprintf('【结论】频域CMA与时域CMA存在差异，需要检查实现\n');
    end
    fprintf('═══════════════════════════════════════════════════════════════\n\n');
    
end

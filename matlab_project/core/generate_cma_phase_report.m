%% generate_cma_phase_report.m: CMA均衡 + 相位噪声补偿综合报告
%
% 功能说明:
%   对TD_TRdata.mat真实数据进行CMA均衡和改进的相位噪声补偿处理
%   生成综合图形化报告，包括:
%   1. 均衡前后星座图对比
%   2. 相位噪声估计与补偿效果
%   3. 误差收敛曲线
%   4. EVM/性能指标统计
%
% 输出:
%   - 综合大图 (整合关键结果)
%   - 可选保存到 reports/ 目录
%
% 作者: DSP课程设计项目
% 日期: 2025年12月

function report_data = generate_cma_phase_report(config)

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
    fprintf('║     CMA均衡 + 相位噪声补偿 综合报告生成器                    ║\n');
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
    
    % 确保输出目录存在
    if config.save_figures && ~exist(config.output_dir, 'dir')
        mkdir(config.output_dir);
    end
    
    % 4QAM星座 (与TD_TRdata数据匹配)
    constellation = [1+1j, 1-1j, -1-1j, -1+1j] / sqrt(2) * sqrt(2);  % R^2 = 2
    
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
    
    S = load(dataFile);
    SigX_in = S.RTdataX(:);
    SigY_in = S.RTdataY(:);
    
    fprintf('  数据长度: %d 样本\n', length(SigX_in));
    
    %% ==================== CMA均衡处理 ====================
    fprintf('执行2x2 MIMO CMA均衡...\n');
    
    % CMA参数
    seg_len = 32;
    tap_len = 33;
    clk_dly = 1;
    par_num = 4;
    step = 2^-8;
    Rx = 2;
    Ry = 2;
    
    seg_num = floor(length(SigX_in) / seg_len);
    clk_num = ceil(seg_num / par_num);
    total_len = seg_len * seg_num;
    
    % 截断数据
    SigX_in = SigX_in(1:total_len);
    SigY_in = SigY_in(1:total_len);
    
    % 前导零填充
    SigX_in_padded = [zeros(tap_len - 1, 1); SigX_in];
    SigY_in_padded = [zeros(tap_len - 1, 1); SigY_in];
    
    % 初始化抽头
    xx = zeros(tap_len, clk_num + clk_dly);
    yy = zeros(tap_len, clk_num + clk_dly);
    xy = zeros(tap_len, clk_num + clk_dly);
    yx = zeros(tap_len, clk_num + clk_dly);
    
    center_tap = (tap_len + 1) / 2;
    xx(center_tap, :) = 1;
    yy(center_tap, :) = 1;
    
    % 输出缓存
    SigX_out = zeros(total_len, 1);
    SigY_out = zeros(total_len, 1);
    errx_vec = zeros(total_len, 1);
    erry_vec = zeros(total_len, 1);
    
    % 抽头累积
    xx_accu = zeros(tap_len, 1);
    yy_accu = zeros(tap_len, 1);
    xy_accu = zeros(tap_len, 1);
    yx_accu = zeros(tap_len, 1);
    
    a = tap_len - 1;
    
    % CMA主循环
    for k = 1:seg_num
        clk_pos = ceil(k / par_num);
        
        xx_curr = xx(:, clk_pos);
        yy_curr = yy(:, clk_pos);
        xy_curr = xy(:, clk_pos);
        yx_curr = yx(:, clk_pos);
        
        for i = 1:seg_len
            current_idx = (k - 1) * seg_len + i;
            padded_idx = current_idx + a;
            input_x = SigX_in_padded(padded_idx:-1:padded_idx - tap_len + 1);
            input_y = SigY_in_padded(padded_idx:-1:padded_idx - tap_len + 1);
            
            SigX_out(current_idx) = xx_curr.' * input_x + xy_curr.' * input_y;
            SigY_out(current_idx) = yx_curr.' * input_x + yy_curr.' * input_y;
            
            errx_vec(current_idx) = Rx - abs(SigX_out(current_idx))^2;
            erry_vec(current_idx) = Ry - abs(SigY_out(current_idx))^2;
            
            xx_accu = xx_accu + step * errx_vec(current_idx) * SigX_out(current_idx) * conj(input_x);
            yy_accu = yy_accu + step * erry_vec(current_idx) * SigY_out(current_idx) * conj(input_y);
            xy_accu = xy_accu + step * errx_vec(current_idx) * SigX_out(current_idx) * conj(input_y);
            yx_accu = yx_accu + step * erry_vec(current_idx) * SigY_out(current_idx) * conj(input_x);
        end
        
        if mod(k, par_num) == 0
            update_clk_pos = clk_pos + clk_dly;
            xx(:, update_clk_pos) = xx(:, update_clk_pos - 1) + xx_accu / par_num;
            yy(:, update_clk_pos) = yy(:, update_clk_pos - 1) + yy_accu / par_num;
            xy(:, update_clk_pos) = xy(:, update_clk_pos - 1) + xy_accu / par_num;
            yx(:, update_clk_pos) = yx(:, update_clk_pos - 1) + yx_accu / par_num;
            
            xx_accu(:) = 0;
            yy_accu(:) = 0;
            xy_accu(:) = 0;
            yx_accu(:) = 0;
        end
        
        if mod(k, 500) == 0
            fprintf('  CMA处理进度: %d / %d (%.1f%%)\n', k, seg_num, 100*k/seg_num);
        end
    end
    
    fprintf('CMA均衡完成.\n');
    
    CMAdataX_raw = SigX_out;
    CMAdataY_raw = SigY_out;
    
    %% ==================== 相位噪声补偿 ====================
    fprintf('执行改进的相位噪声补偿...\n');
    
    % 使用改进的相位噪声估计器
    pn_config = struct();
    pn_config.mod_order = 4;
    pn_config.win_half = 5;
    pn_config.constellation = constellation;
    pn_config.ambiguity_method = 'unwrap';
    pn_config.phase_unwrap = true;
    
    [phase_est_X, CMAdataX_corr, diagX] = phase_noise_estimator(CMAdataX_raw, pn_config);
    [phase_est_Y, CMAdataY_corr, diagY] = phase_noise_estimator(CMAdataY_raw, pn_config);
    
    fprintf('  X极化: EVM %.2f%% -> %.2f%%\n', diagX.evm_before, diagX.evm_after);
    fprintf('  Y极化: EVM %.2f%% -> %.2f%%\n', diagY.evm_before, diagY.evm_after);
    
    %% ==================== 综合大图1: 处理流程对比 ====================
    fprintf('\n生成综合报告图1: 处理流程对比...\n');
    
    figure('Name', 'CMA+相位补偿综合报告', 'Position', [50, 50, 1600, 900]);
    
    % 取末尾稳定部分用于显示
    num_show = min(50000, total_len - 1);
    idx_show = (total_len - num_show):total_len;
    
    % 第一行: X极化处理流程
    subplot(3, 4, 1);
    plot(real(SigX_in(idx_show)), imag(SigX_in(idx_show)), '.r', 'MarkerSize', 1);
    title('X极化 - 均衡前');
    xlabel('实部'); ylabel('虚部');
    axis equal; grid on;
    xlim([-2.5 2.5]); ylim([-2.5 2.5]);
    
    subplot(3, 4, 2);
    plot(real(CMAdataX_raw(idx_show)), imag(CMAdataX_raw(idx_show)), '.', 'Color', [1 0.5 0], 'MarkerSize', 1);
    title(sprintf('X极化 - CMA均衡后\nEVM=%.1f%%', diagX.evm_before));
    xlabel('实部'); ylabel('虚部');
    axis equal; grid on;
    xlim([-2.5 2.5]); ylim([-2.5 2.5]);
    
    subplot(3, 4, 3);
    plot(real(CMAdataX_corr(idx_show)), imag(CMAdataX_corr(idx_show)), '.g', 'MarkerSize', 1);
    hold on;
    plot(real(constellation), imag(constellation), 'ko', 'MarkerSize', 12, 'LineWidth', 2);
    title(sprintf('X极化 - 相位补偿后\nEVM=%.1f%%', diagX.evm_after));
    xlabel('实部'); ylabel('虚部');
    axis equal; grid on;
    xlim([-2.5 2.5]); ylim([-2.5 2.5]);
    
    % X极化相位轨迹 (后100点)
    subplot(3, 4, 4);
    num_phase_show = 100;
    idx_phase_show = (total_len - num_phase_show + 1):total_len;
    phase_X_plot = mod(phase_est_X(idx_phase_show), 2*pi) * 180/pi;
    plot(1:num_phase_show, phase_X_plot, 'b-o', 'LineWidth', 1, 'MarkerSize', 3);
    xlabel('符号索引 (后100点)'); ylabel('相位 (度)');
    title('X极化 - 估计相位轨迹');
    grid on;
    ylim([0 360]);
    
    % 第二行: Y极化处理流程
    subplot(3, 4, 5);
    plot(real(SigY_in(idx_show)), imag(SigY_in(idx_show)), '.b', 'MarkerSize', 1);
    title('Y极化 - 均衡前');
    xlabel('实部'); ylabel('虚部');
    axis equal; grid on;
    xlim([-2.5 2.5]); ylim([-2.5 2.5]);
    
    subplot(3, 4, 6);
    plot(real(CMAdataY_raw(idx_show)), imag(CMAdataY_raw(idx_show)), '.', 'Color', [1 0.5 0], 'MarkerSize', 1);
    title(sprintf('Y极化 - CMA均衡后\nEVM=%.1f%%', diagY.evm_before));
    xlabel('实部'); ylabel('虚部');
    axis equal; grid on;
    xlim([-2.5 2.5]); ylim([-2.5 2.5]);
    
    subplot(3, 4, 7);
    plot(real(CMAdataY_corr(idx_show)), imag(CMAdataY_corr(idx_show)), '.c', 'MarkerSize', 1);
    hold on;
    plot(real(constellation), imag(constellation), 'ko', 'MarkerSize', 12, 'LineWidth', 2);
    title(sprintf('Y极化 - 相位补偿后\nEVM=%.1f%%', diagY.evm_after));
    xlabel('实部'); ylabel('虚部');
    axis equal; grid on;
    xlim([-2.5 2.5]); ylim([-2.5 2.5]);
    
    % Y极化相位轨迹 (后100点)
    subplot(3, 4, 8);
    phase_Y_plot = mod(phase_est_Y(idx_phase_show), 2*pi) * 180/pi;
    plot(1:num_phase_show, phase_Y_plot, 'r-o', 'LineWidth', 1, 'MarkerSize', 3);
    xlabel('符号索引 (后100点)'); ylabel('相位 (度)');
    title('Y极化 - 估计相位轨迹');
    grid on;
    ylim([0 360]);
    
    % 第三行: 误差收敛和性能统计
    subplot(3, 4, 9);
    errx_smooth = movmean(abs(errx_vec), 1000);
    erry_smooth = movmean(abs(erry_vec), 1000);
    plot(errx_smooth, 'r-', 'LineWidth', 1);
    hold on;
    plot(erry_smooth, 'b-', 'LineWidth', 1);
    xlabel('符号索引'); ylabel('|误差|');
    title('CMA误差收敛曲线');
    legend('X极化', 'Y极化', 'Location', 'northeast');
    grid on;
    
    subplot(3, 4, 10);
    % 误差直方图
    histogram(abs(errx_vec(end-50000:end)), 50, 'FaceColor', 'r', 'FaceAlpha', 0.5);
    hold on;
    histogram(abs(erry_vec(end-50000:end)), 50, 'FaceColor', 'b', 'FaceAlpha', 0.5);
    xlabel('|误差|'); ylabel('计数');
    title('稳态误差分布');
    legend('X极化', 'Y极化');
    grid on;
    
    % 性能表格
    subplot(3, 4, 11);
    axis off;
    text(0.1, 0.9, '【性能指标摘要】', 'FontSize', 12, 'FontWeight', 'bold');
    text(0.1, 0.75, sprintf('数据长度: %d 符号', total_len), 'FontSize', 10);
    text(0.1, 0.60, sprintf('CMA抽头数: %d', tap_len), 'FontSize', 10);
    text(0.1, 0.45, sprintf('X极化 EVM改善: %.1f%% → %.1f%%', diagX.evm_before, diagX.evm_after), 'FontSize', 10);
    text(0.1, 0.30, sprintf('Y极化 EVM改善: %.1f%% → %.1f%%', diagY.evm_before, diagY.evm_after), 'FontSize', 10);
    text(0.1, 0.15, sprintf('相位噪声估计窗口: %d', 2*pn_config.win_half+1), 'FontSize', 10);
    text(0.1, 0.0, sprintf('处理时间: %s', datestr(now, 'yyyy-mm-dd HH:MM:SS')), 'FontSize', 10);
    
    % 相位补偿方法说明
    subplot(3, 4, 12);
    axis off;
    text(0.1, 0.9, '【算法说明】', 'FontSize', 12, 'FontWeight', 'bold');
    text(0.1, 0.75, '1. 2×2 MIMO CMA盲均衡', 'FontSize', 10);
    text(0.1, 0.60, '2. 四次方V&V相位估计', 'FontSize', 10);
    text(0.1, 0.45, '3. 星座固有相位补偿', 'FontSize', 10);
    text(0.1, 0.30, '4. 相位模糊自动解决', 'FontSize', 10);
    text(0.1, 0.15, '调制: 4QAM/QPSK', 'FontSize', 10);
    text(0.1, 0.0, sprintf('模糊解决方法: %s', pn_config.ambiguity_method), 'FontSize', 10);
    
    sgtitle('【综合报告】CMA均衡 + 相位噪声补偿 - TD\_TRdata处理结果', 'FontSize', 14, 'FontWeight', 'bold');
    
    if config.save_figures
        savefig(fullfile(config.output_dir, 'cma_phase_report_main.fig'));
        print(fullfile(config.output_dir, 'cma_phase_report_main.png'), '-dpng', '-r150');
        fprintf('  保存到: %s\n', fullfile(config.output_dir, 'cma_phase_report_main.png'));
    end
    
    %% ==================== 综合大图2: 相位跟踪详细分析 ====================
    fprintf('生成综合报告图2: 相位跟踪详细分析...\n');
    
    figure('Name', '相位跟踪详细分析', 'Position', [100, 100, 1400, 800]);
    
    % 取后100点数据进行详细分析，便于看清相位跟随效果
    num_detail = 100;
    idx_detail = (total_len - num_detail):total_len;
    
    % 相位轨迹对比
    subplot(2, 3, 1);
    phase_X_detail = mod(phase_est_X(idx_detail), 2*pi) * 180/pi;
    phase_Y_detail = mod(phase_est_Y(idx_detail), 2*pi) * 180/pi;
    plot(1:num_detail+1, phase_X_detail, 'b-', 'LineWidth', 1);
    hold on;
    plot(1:num_detail+1, phase_Y_detail, 'r-', 'LineWidth', 1);
    xlabel('相对符号索引');
    ylabel('相位 (度)');
    title('X/Y极化相位轨迹对比 (0°~360°)');
    legend('X极化', 'Y极化', 'Location', 'best');
    grid on;
    ylim([0 360]);
    
    % 相位差
    subplot(2, 3, 2);
    phase_diff = phase_X_detail - phase_Y_detail;
    phase_diff = mod(phase_diff + 180, 360) - 180;  % 归一化到 [-180, 180]
    plot(1:num_detail+1, phase_diff, 'k-', 'LineWidth', 0.5);
    xlabel('相对符号索引');
    ylabel('相位差 (度)');
    title(sprintf('X-Y极化相位差 (均值=%.1f°, std=%.1f°)', mean(phase_diff), std(phase_diff)));
    grid on;
    ylim([-180 180]);
    
    % 相位变化率
    subplot(2, 3, 3);
    phase_rate_X = diff(unwrap(phase_est_X(idx_detail))) * 180/pi;
    phase_rate_Y = diff(unwrap(phase_est_Y(idx_detail))) * 180/pi;
    histogram(phase_rate_X, 100, 'FaceColor', 'b', 'FaceAlpha', 0.5, 'Normalization', 'probability');
    hold on;
    histogram(phase_rate_Y, 100, 'FaceColor', 'r', 'FaceAlpha', 0.5, 'Normalization', 'probability');
    xlabel('相位变化率 (度/符号)');
    ylabel('概率');
    title(sprintf('相位变化率分布\nX: σ=%.2f°, Y: σ=%.2f°', std(phase_rate_X), std(phase_rate_Y)));
    legend('X极化', 'Y极化');
    grid on;
    xlim([-10 10]);
    
    % 补偿前后星座图 (使用全部点)
    subplot(2, 3, 4);
    scatter(real(CMAdataX_raw), imag(CMAdataX_raw), 1, 'r', 'filled', 'MarkerFaceAlpha', 0.1);
    hold on;
    scatter(real(CMAdataX_corr), imag(CMAdataX_corr), 1, 'b', 'filled', 'MarkerFaceAlpha', 0.1);
    plot(real(constellation), imag(constellation), 'ko', 'MarkerSize', 15, 'LineWidth', 2);
    xlabel('实部'); ylabel('虚部');
    title('X极化: 相位补偿前(红) vs 补偿后(蓝)');
    axis equal; grid on;
    xlim([-2.5 2.5]); ylim([-2.5 2.5]);
    legend('补偿前', '补偿后', '理想星座', 'Location', 'best');
    
    subplot(2, 3, 5);
    scatter(real(CMAdataY_raw), imag(CMAdataY_raw), 1, 'r', 'filled', 'MarkerFaceAlpha', 0.1);
    hold on;
    scatter(real(CMAdataY_corr), imag(CMAdataY_corr), 1, 'b', 'filled', 'MarkerFaceAlpha', 0.1);
    plot(real(constellation), imag(constellation), 'ko', 'MarkerSize', 15, 'LineWidth', 2);
    xlabel('实部'); ylabel('虚部');
    title('Y极化: 相位补偿前(红) vs 补偿后(蓝)');
    axis equal; grid on;
    xlim([-2.5 2.5]); ylim([-2.5 2.5]);
    legend('补偿前', '补偿后', '理想星座', 'Location', 'best');
    
    % 综合性能柱状图
    subplot(2, 3, 6);
    categories = {'X极化补偿前', 'X极化补偿后', 'Y极化补偿前', 'Y极化补偿后'};
    evm_values = [diagX.evm_before, diagX.evm_after, diagY.evm_before, diagY.evm_after];
    bar_colors = [1 0.3 0.3; 0.3 0.7 0.3; 1 0.3 0.3; 0.3 0.7 0.3];
    b = bar(evm_values);
    b.FaceColor = 'flat';
    b.CData = bar_colors;
    set(gca, 'XTickLabel', categories);
    ylabel('EVM (%)');
    title('EVM性能对比');
    grid on;
    
    % 添加数值标签
    for i = 1:length(evm_values)
        text(i, evm_values(i) + 1, sprintf('%.1f%%', evm_values(i)), ...
             'HorizontalAlignment', 'center', 'FontSize', 10);
    end
    
    sgtitle('【详细分析】相位跟踪与补偿效果', 'FontSize', 14, 'FontWeight', 'bold');
    
    if config.save_figures
        savefig(fullfile(config.output_dir, 'cma_phase_report_detail.fig'));
        print(fullfile(config.output_dir, 'cma_phase_report_detail.png'), '-dpng', '-r150');
        fprintf('  保存到: %s\n', fullfile(config.output_dir, 'cma_phase_report_detail.png'));
    end
    
    %% ==================== 保存报告数据 ====================
    report_data.SigX_in = SigX_in;
    report_data.SigY_in = SigY_in;
    report_data.CMAdataX_raw = CMAdataX_raw;
    report_data.CMAdataY_raw = CMAdataY_raw;
    report_data.CMAdataX_corr = CMAdataX_corr;
    report_data.CMAdataY_corr = CMAdataY_corr;
    report_data.phase_est_X = phase_est_X;
    report_data.phase_est_Y = phase_est_Y;
    report_data.errx_vec = errx_vec;
    report_data.erry_vec = erry_vec;
    report_data.diagX = diagX;
    report_data.diagY = diagY;
    report_data.constellation = constellation;
    
    if config.save_figures
        save(fullfile(config.output_dir, 'cma_phase_report_data.mat'), 'report_data');
        fprintf('  数据保存到: %s\n', fullfile(config.output_dir, 'cma_phase_report_data.mat'));
    end
    
    fprintf('\n报告生成完成!\n');
    fprintf('═══════════════════════════════════════════════════════════════\n');
    fprintf('  X极化 EVM: %.2f%% → %.2f%% (改善 %.1f%%)\n', ...
            diagX.evm_before, diagX.evm_after, diagX.evm_before - diagX.evm_after);
    fprintf('  Y极化 EVM: %.2f%% → %.2f%% (改善 %.1f%%)\n', ...
            diagY.evm_before, diagY.evm_after, diagY.evm_before - diagY.evm_after);
    fprintf('═══════════════════════════════════════════════════════════════\n\n');
    
end

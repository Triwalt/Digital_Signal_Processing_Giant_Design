%% generate_cma_parameter_study.m: CMA参数敏感性分析报告
%
% 功能说明:
%   评估CMA均衡器关键参数对均衡效果的影响:
%   1. 步长 (step) - 影响收敛速度和稳态误差
%   2. 抽头长度 (tap_len) - 影响均衡能力和计算复杂度
%   3. 分段长度 (seg_len) - 影响延迟和并行效率
%   4. 并行段数 (par_num) - 影响抽头更新频率
%
% 输出:
%   - 综合报告图 (PNG)
%   - 参数扫描数据 (MAT)
%
% 作者: DSP课程设计项目
% 日期: 2025年12月

function report_data = generate_cma_parameter_study(config)

    clc;
    close all;
    
    %% 路径设置
    scriptDir = fileparts(mfilename('fullpath'));
    if isempty(scriptDir)
        scriptDir = pwd;
    end
    
    % 添加路径
    coreDir = fullfile(scriptDir);
    hwDir = fullfile(scriptDir, '..', 'homework1');
    addpath(coreDir);
    addpath(hwDir);
    
    fprintf('╔══════════════════════════════════════════════════════════════╗\n');
    fprintf('║     CMA均衡器参数敏感性分析报告生成器                        ║\n');
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
    
    %% 加载测试数据
    fprintf('加载TD_TRdata.mat数据...\n');
    dataFile = fullfile(hwDir, 'TD_TRdata.mat');
    S = load(dataFile);
    SigX_in = S.RTdataX(:);
    SigY_in = S.RTdataY(:);
    fprintf('  数据长度: %d 样本\n\n', length(SigX_in));
    
    % 报告数据
    report_data = struct();
    report_data.config = config;
    report_data.timestamp = datetime('now');
    
    %% ==================== 基准测试 ====================
    fprintf('[0/4] 运行基准配置...\n');
    
    base_params = struct();
    base_params.step = 2^-8;
    base_params.tap_len = 33;
    base_params.seg_len = 32;
    base_params.par_num = 4;
    base_params.use_builtin_fft = true;
    
    base_result = cma_equalizer_core(SigX_in, SigY_in, base_params);
    fprintf('  基准EVM: %.2f%% -> %.2f%%\n', base_result.evm_before, base_result.evm_after);
    fprintf('  基准收敛点: %d 样本\n', base_result.convergence_idx);
    fprintf('  基准执行时间: %.3f 秒\n\n', base_result.exec_time);
    
    report_data.base_params = base_params;
    report_data.base_result = base_result;
    
    %% ==================== 1. 步长扫描 ====================
    fprintf('[1/4] 步长扫描...\n');
    
    step_values = [2^-12, 2^-11, 2^-10, 2^-9, 2^-8, 2^-7, 2^-6, 2^-5, 2^-4];
    step_labels = {'2^{-12}', '2^{-11}', '2^{-10}', '2^{-9}', '2^{-8}', '2^{-7}', '2^{-6}', '2^{-5}', '2^{-4}'};
    n_step = length(step_values);
    
    step_study = struct();
    step_study.values = step_values;
    step_study.labels = step_labels;
    step_study.evm_after = zeros(n_step, 1);
    step_study.convergence_idx = zeros(n_step, 1);
    step_study.exec_time = zeros(n_step, 1);
    step_study.steady_state_err = zeros(n_step, 1);
    step_study.results = cell(n_step, 1);
    
    for i = 1:n_step
        params = base_params;
        params.step = step_values(i);
        
        result = cma_equalizer_core(SigX_in, SigY_in, params);
        
        step_study.evm_after(i) = result.evm_after;
        step_study.convergence_idx(i) = result.convergence_idx;
        step_study.exec_time(i) = result.exec_time;
        step_study.steady_state_err(i) = mean(abs(result.errx_vec(end-10000:end)));
        step_study.results{i} = result;
        
        fprintf('  step=%s: EVM=%.2f%%, 收敛=%d, 稳态误差=%.4f\n', ...
                step_labels{i}, result.evm_after, result.convergence_idx, step_study.steady_state_err(i));
    end
    
    report_data.step_study = step_study;
    
    %% ==================== 2. 抽头长度扫描 ====================
    fprintf('\n[2/4] 抽头长度扫描...\n');
    
    tap_values = [9, 17, 25, 33, 49, 65, 97, 129];
    n_tap = length(tap_values);
    
    tap_study = struct();
    tap_study.values = tap_values;
    tap_study.evm_after = zeros(n_tap, 1);
    tap_study.convergence_idx = zeros(n_tap, 1);
    tap_study.exec_time = zeros(n_tap, 1);
    tap_study.results = cell(n_tap, 1);
    
    for i = 1:n_tap
        params = base_params;
        params.tap_len = tap_values(i);
        
        result = cma_equalizer_core(SigX_in, SigY_in, params);
        
        tap_study.evm_after(i) = result.evm_after;
        tap_study.convergence_idx(i) = result.convergence_idx;
        tap_study.exec_time(i) = result.exec_time;
        tap_study.results{i} = result;
        
        fprintf('  tap_len=%d: EVM=%.2f%%, 收敛=%d, 时间=%.3fs\n', ...
                tap_values(i), result.evm_after, result.convergence_idx, result.exec_time);
    end
    
    report_data.tap_study = tap_study;
    
    %% ==================== 3. 分段长度扫描 ====================
    fprintf('\n[3/4] 分段长度扫描...\n');
    
    seg_values = [8, 16, 32, 64, 128, 256];
    n_seg = length(seg_values);
    
    seg_study = struct();
    seg_study.values = seg_values;
    seg_study.evm_after = zeros(n_seg, 1);
    seg_study.convergence_idx = zeros(n_seg, 1);
    seg_study.exec_time = zeros(n_seg, 1);
    seg_study.results = cell(n_seg, 1);
    
    for i = 1:n_seg
        params = base_params;
        params.seg_len = seg_values(i);
        
        result = cma_equalizer_core(SigX_in, SigY_in, params);
        
        seg_study.evm_after(i) = result.evm_after;
        seg_study.convergence_idx(i) = result.convergence_idx;
        seg_study.exec_time(i) = result.exec_time;
        seg_study.results{i} = result;
        
        fprintf('  seg_len=%d: EVM=%.2f%%, 收敛=%d, 时间=%.3fs\n', ...
                seg_values(i), result.evm_after, result.convergence_idx, result.exec_time);
    end
    
    report_data.seg_study = seg_study;
    
    %% ==================== 4. 并行段数扫描 ====================
    fprintf('\n[4/4] 并行段数扫描...\n');
    
    par_values = [1, 2, 4, 8, 16, 32];
    n_par = length(par_values);
    
    par_study = struct();
    par_study.values = par_values;
    par_study.evm_after = zeros(n_par, 1);
    par_study.convergence_idx = zeros(n_par, 1);
    par_study.exec_time = zeros(n_par, 1);
    par_study.results = cell(n_par, 1);
    
    for i = 1:n_par
        params = base_params;
        params.par_num = par_values(i);
        
        result = cma_equalizer_core(SigX_in, SigY_in, params);
        
        par_study.evm_after(i) = result.evm_after;
        par_study.convergence_idx(i) = result.convergence_idx;
        par_study.exec_time(i) = result.exec_time;
        par_study.results{i} = result;
        
        fprintf('  par_num=%d: EVM=%.2f%%, 收敛=%d, 时间=%.3fs\n', ...
                par_values(i), result.evm_after, result.convergence_idx, result.exec_time);
    end
    
    report_data.par_study = par_study;
    
    %% ==================== 生成综合报告图1: 参数扫描结果 ====================
    fprintf('\n生成综合报告图1...\n');
    
    figure('Name', 'CMA参数敏感性分析', 'Position', [50, 50, 1600, 900]);
    
    % 步长分析
    subplot(3, 4, 1);
    semilogx(step_values, step_study.evm_after, 'b-o', 'LineWidth', 2, 'MarkerSize', 8);
    xlabel('步长'); ylabel('EVM (%)');
    title('步长 vs EVM');
    grid on;
    set(gca, 'XTick', step_values(1:2:end));
    
    subplot(3, 4, 2);
    semilogx(step_values, step_study.convergence_idx / 1000, 'r-s', 'LineWidth', 2, 'MarkerSize', 8);
    xlabel('步长'); ylabel('收敛点 (×1000)');
    title('步长 vs 收敛速度');
    grid on;
    set(gca, 'XTick', step_values(1:2:end));
    
    subplot(3, 4, 3);
    semilogx(step_values, step_study.steady_state_err, 'g-^', 'LineWidth', 2, 'MarkerSize', 8);
    xlabel('步长'); ylabel('稳态误差');
    title('步长 vs 稳态误差');
    grid on;
    set(gca, 'XTick', step_values(1:2:end));
    
    % 步长收敛曲线对比
    subplot(3, 4, 4);
    colors = lines(n_step);
    legend_entries = {};
    for i = [1, 3, 5, 7, 9]  % 选择部分显示
        if i <= n_step
            err_smooth = movmean(abs(step_study.results{i}.errx_vec), 1000);
            plot(err_smooth(1:1000:end), 'Color', colors(i,:), 'LineWidth', 1.5);
            hold on;
            legend_entries{end+1} = step_labels{i};
        end
    end
    xlabel('样本索引 (×1000)'); ylabel('|误差|');
    title('不同步长收敛曲线');
    legend(legend_entries, 'Location', 'northeast');
    grid on;
    
    % 抽头长度分析
    subplot(3, 4, 5);
    plot(tap_values, tap_study.evm_after, 'b-o', 'LineWidth', 2, 'MarkerSize', 8);
    xlabel('抽头长度'); ylabel('EVM (%)');
    title('抽头长度 vs EVM');
    grid on;
    
    subplot(3, 4, 6);
    plot(tap_values, tap_study.convergence_idx / 1000, 'r-s', 'LineWidth', 2, 'MarkerSize', 8);
    xlabel('抽头长度'); ylabel('收敛点 (×1000)');
    title('抽头长度 vs 收敛速度');
    grid on;
    
    subplot(3, 4, 7);
    plot(tap_values, tap_study.exec_time, 'g-^', 'LineWidth', 2, 'MarkerSize', 8);
    xlabel('抽头长度'); ylabel('执行时间 (s)');
    title('抽头长度 vs 计算时间');
    grid on;
    
    % 抽头长度星座图对比
    subplot(3, 4, 8);
    idx_tap_small = 1;  % 最小抽头
    idx_tap_large = n_tap;  % 最大抽头
    plot(tap_study.results{idx_tap_small}.SigX_out(end-5000:end), '.r', 'MarkerSize', 2);
    hold on;
    plot(tap_study.results{idx_tap_large}.SigX_out(end-5000:end), '.b', 'MarkerSize', 2);
    title(sprintf('抽头长度对比: %d vs %d', tap_values(idx_tap_small), tap_values(idx_tap_large)));
    xlabel('实部'); ylabel('虚部');
    axis equal; grid on;
    xlim([-2.5 2.5]); ylim([-2.5 2.5]);
    legend(sprintf('tap=%d', tap_values(idx_tap_small)), sprintf('tap=%d', tap_values(idx_tap_large)));
    
    % 分段长度分析
    subplot(3, 4, 9);
    semilogx(seg_values, seg_study.evm_after, 'b-o', 'LineWidth', 2, 'MarkerSize', 8);
    xlabel('分段长度'); ylabel('EVM (%)');
    title('分段长度 vs EVM');
    grid on;
    set(gca, 'XTick', seg_values);
    
    subplot(3, 4, 10);
    semilogx(seg_values, seg_study.exec_time, 'g-^', 'LineWidth', 2, 'MarkerSize', 8);
    xlabel('分段长度'); ylabel('执行时间 (s)');
    title('分段长度 vs 计算时间');
    grid on;
    set(gca, 'XTick', seg_values);
    
    % 并行段数分析
    subplot(3, 4, 11);
    semilogx(par_values, par_study.evm_after, 'b-o', 'LineWidth', 2, 'MarkerSize', 8);
    xlabel('并行段数'); ylabel('EVM (%)');
    title('并行段数 vs EVM');
    grid on;
    set(gca, 'XTick', par_values);
    
    subplot(3, 4, 12);
    semilogx(par_values, par_study.convergence_idx / 1000, 'r-s', 'LineWidth', 2, 'MarkerSize', 8);
    xlabel('并行段数'); ylabel('收敛点 (×1000)');
    title('并行段数 vs 收敛速度');
    grid on;
    set(gca, 'XTick', par_values);
    
    sgtitle('【综合报告】CMA均衡器参数敏感性分析', 'FontSize', 14, 'FontWeight', 'bold');
    
    if config.save_figures
        savefig(fullfile(config.output_dir, 'cma_parameter_study.fig'));
        print(fullfile(config.output_dir, 'cma_parameter_study.png'), '-dpng', '-r150');
        fprintf('  保存到: %s\n', fullfile(config.output_dir, 'cma_parameter_study.png'));
    end
    
    %% ==================== 生成综合报告图2: 最优参数建议 ====================
    fprintf('生成综合报告图2...\n');
    
    figure('Name', 'CMA参数优化建议', 'Position', [100, 100, 1400, 800]);
    
    % 找出最优参数
    [~, best_step_idx] = min(step_study.evm_after);
    [~, best_tap_idx] = min(tap_study.evm_after);
    [~, best_seg_idx] = min(seg_study.evm_after);
    [~, best_par_idx] = min(par_study.evm_after);
    
    % 综合评分 (EVM权重0.5, 收敛速度0.3, 计算时间0.2)
    step_score = normalize_score(step_study.evm_after, 0.5) + ...
                 normalize_score(step_study.convergence_idx, 0.3) + ...
                 normalize_score(step_study.exec_time, 0.2);
    [~, best_step_overall] = min(step_score);
    
    tap_score = normalize_score(tap_study.evm_after, 0.5) + ...
                normalize_score(tap_study.convergence_idx, 0.3) + ...
                normalize_score(tap_study.exec_time, 0.2);
    [~, best_tap_overall] = min(tap_score);
    
    % 步长权衡图
    subplot(2, 3, 1);
    yyaxis left;
    semilogx(step_values, step_study.evm_after, 'b-o', 'LineWidth', 2, 'MarkerSize', 8);
    ylabel('EVM (%)');
    yyaxis right;
    semilogx(step_values, step_study.convergence_idx / 1000, 'r-s', 'LineWidth', 2, 'MarkerSize', 8);
    ylabel('收敛点 (×1000)');
    xlabel('步长');
    title(sprintf('步长权衡 (推荐: %s)', step_labels{best_step_overall}));
    grid on;
    xline(step_values(best_step_overall), 'g--', 'LineWidth', 2);
    
    % 抽头长度权衡图
    subplot(2, 3, 2);
    yyaxis left;
    plot(tap_values, tap_study.evm_after, 'b-o', 'LineWidth', 2, 'MarkerSize', 8);
    ylabel('EVM (%)');
    yyaxis right;
    plot(tap_values, tap_study.exec_time, 'r-s', 'LineWidth', 2, 'MarkerSize', 8);
    ylabel('执行时间 (s)');
    xlabel('抽头长度');
    title(sprintf('抽头长度权衡 (推荐: %d)', tap_values(best_tap_overall)));
    grid on;
    xline(tap_values(best_tap_overall), 'g--', 'LineWidth', 2);
    
    % 星座图对比: 基准 vs 最优
    subplot(2, 3, 3);
    plot(base_result.SigX_out(end-10000:end), '.r', 'MarkerSize', 1);
    title(sprintf('基准配置星座图\nEVM=%.2f%%', base_result.evm_after));
    xlabel('实部'); ylabel('虚部');
    axis equal; grid on;
    xlim([-2.5 2.5]); ylim([-2.5 2.5]);
    
    % 最优配置测试
    subplot(2, 3, 4);
    optimal_params = base_params;
    optimal_params.step = step_values(best_step_overall);
    optimal_params.tap_len = tap_values(best_tap_overall);
    optimal_result = cma_equalizer_core(SigX_in, SigY_in, optimal_params);
    plot(optimal_result.SigX_out(end-10000:end), '.b', 'MarkerSize', 1);
    title(sprintf('优化配置星座图\nEVM=%.2f%%', optimal_result.evm_after));
    xlabel('实部'); ylabel('虚部');
    axis equal; grid on;
    xlim([-2.5 2.5]); ylim([-2.5 2.5]);
    
    report_data.optimal_params = optimal_params;
    report_data.optimal_result = optimal_result;
    
    % 参数影响热力图
    subplot(2, 3, 5);
    % 创建步长-抽头二维扫描的简化版本
    step_subset = [2^-10, 2^-8, 2^-6];
    tap_subset = [17, 33, 65];
    evm_matrix = zeros(length(step_subset), length(tap_subset));
    
    for i = 1:length(step_subset)
        for j = 1:length(tap_subset)
            params = base_params;
            params.step = step_subset(i);
            params.tap_len = tap_subset(j);
            result = cma_equalizer_core(SigX_in, SigY_in, params);
            evm_matrix(i, j) = result.evm_after;
        end
    end
    
    imagesc(evm_matrix);
    colorbar;
    colormap(flipud(hot));
    set(gca, 'XTick', 1:length(tap_subset), 'XTickLabel', tap_subset);
    set(gca, 'YTick', 1:length(step_subset), 'YTickLabel', {'2^{-10}', '2^{-8}', '2^{-6}'});
    xlabel('抽头长度'); ylabel('步长');
    title('步长-抽头长度 EVM热力图 (%)');
    
    % 参数建议表格
    subplot(2, 3, 6);
    axis off;
    text(0.05, 0.95, '【参数优化建议】', 'FontSize', 14, 'FontWeight', 'bold');
    text(0.05, 0.82, '─────────────────────────────────', 'FontSize', 10);
    text(0.05, 0.72, sprintf('基准配置:'), 'FontSize', 11, 'FontWeight', 'bold');
    text(0.05, 0.62, sprintf('  步长: 2^{-8}, 抽头: 33, EVM: %.2f%%', base_result.evm_after), 'FontSize', 10);
    text(0.05, 0.50, sprintf('推荐配置:'), 'FontSize', 11, 'FontWeight', 'bold', 'Color', 'b');
    text(0.05, 0.40, sprintf('  步长: %s', step_labels{best_step_overall}), 'FontSize', 10);
    text(0.05, 0.30, sprintf('  抽头长度: %d', tap_values(best_tap_overall)), 'FontSize', 10);
    text(0.05, 0.20, sprintf('  优化后EVM: %.2f%%', optimal_result.evm_after), 'FontSize', 10);
    text(0.05, 0.08, '─────────────────────────────────', 'FontSize', 10);
    
    % 关键发现
    if step_study.evm_after(1) > step_study.evm_after(end)
        finding1 = '步长↑ → EVM↓ (但稳态误差↑)';
    else
        finding1 = '步长↓ → EVM↓ (收敛更慢)';
    end
    if tap_study.evm_after(1) > tap_study.evm_after(end)
        finding2 = '抽头↑ → EVM↓ (计算量↑)';
    else
        finding2 = '抽头影响不显著';
    end
    text(0.05, -0.05, sprintf('发现1: %s', finding1), 'FontSize', 9);
    text(0.05, -0.15, sprintf('发现2: %s', finding2), 'FontSize', 9);
    
    sgtitle('【优化报告】CMA均衡器最优参数建议', 'FontSize', 14, 'FontWeight', 'bold');
    
    if config.save_figures
        savefig(fullfile(config.output_dir, 'cma_parameter_optimal.fig'));
        print(fullfile(config.output_dir, 'cma_parameter_optimal.png'), '-dpng', '-r150');
        fprintf('  保存到: %s\n', fullfile(config.output_dir, 'cma_parameter_optimal.png'));
    end
    
    %% ==================== 保存报告数据 ====================
    if config.save_figures
        save(fullfile(config.output_dir, 'cma_parameter_study_data.mat'), 'report_data');
        fprintf('  数据保存到: %s\n', fullfile(config.output_dir, 'cma_parameter_study_data.mat'));
    end
    
    %% ==================== 打印最终报告 ====================
    fprintf('\n');
    fprintf('═══════════════════════════════════════════════════════════════\n');
    fprintf('                CMA参数敏感性分析报告                          \n');
    fprintf('═══════════════════════════════════════════════════════════════\n');
    fprintf('【步长分析】 (基准: 2^-8)\n');
    fprintf('  最小EVM: %.2f%% @ step=%s\n', min(step_study.evm_after), step_labels{best_step_idx});
    fprintf('  最快收敛: %d @ step=%s\n', min(step_study.convergence_idx), step_labels{find(step_study.convergence_idx == min(step_study.convergence_idx), 1)});
    fprintf('  建议: %s (综合评分最优)\n', step_labels{best_step_overall});
    fprintf('───────────────────────────────────────────────────────────────\n');
    fprintf('【抽头长度分析】 (基准: 33)\n');
    fprintf('  最小EVM: %.2f%% @ tap=%d\n', min(tap_study.evm_after), tap_values(best_tap_idx));
    fprintf('  最快执行: %.3fs @ tap=%d\n', min(tap_study.exec_time), tap_values(find(tap_study.exec_time == min(tap_study.exec_time), 1)));
    fprintf('  建议: %d (综合评分最优)\n', tap_values(best_tap_overall));
    fprintf('───────────────────────────────────────────────────────────────\n');
    fprintf('【分段长度分析】 (基准: 32)\n');
    fprintf('  EVM范围: %.2f%% ~ %.2f%%\n', min(seg_study.evm_after), max(seg_study.evm_after));
    fprintf('  执行时间范围: %.3fs ~ %.3fs\n', min(seg_study.exec_time), max(seg_study.exec_time));
    fprintf('───────────────────────────────────────────────────────────────\n');
    fprintf('【并行段数分析】 (基准: 4)\n');
    fprintf('  EVM范围: %.2f%% ~ %.2f%%\n', min(par_study.evm_after), max(par_study.evm_after));
    fprintf('───────────────────────────────────────────────────────────────\n');
    fprintf('【最优配置】\n');
    fprintf('  步长: %s, 抽头: %d\n', step_labels{best_step_overall}, tap_values(best_tap_overall));
    fprintf('  EVM: %.2f%% (改善 %.2f%%)\n', optimal_result.evm_after, base_result.evm_after - optimal_result.evm_after);
    fprintf('═══════════════════════════════════════════════════════════════\n\n');
    
end

%% 辅助函数: 归一化评分
function score = normalize_score(values, weight)
    % 将值归一化到 [0, 1] 范围，值越小得分越低
    min_val = min(values);
    max_val = max(values);
    if max_val == min_val
        score = zeros(size(values));
    else
        score = (values - min_val) / (max_val - min_val) * weight;
    end
end

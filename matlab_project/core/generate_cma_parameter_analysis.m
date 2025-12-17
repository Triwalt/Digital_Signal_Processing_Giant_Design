%% generate_cma_parameter_analysis.m: CMA参数性能深度分析报告
%
% 功能说明:
%   基于参数扫描结果，深入分析各参数对CMA性能的影响机理
%   包括理论分析、趋势解释和参数调优建议
%
% 作者: DSP课程设计项目
% 日期: 2025年12月

function generate_cma_parameter_analysis()

    clc;
    close all;
    
    %% 路径设置
    scriptDir = fileparts(mfilename('fullpath'));
    if isempty(scriptDir)
        scriptDir = pwd;
    end
    
    reportDir = fullfile(scriptDir, '..', 'reports');
    
    fprintf('╔══════════════════════════════════════════════════════════════╗\n');
    fprintf('║     CMA均衡器参数性能深度分析报告                            ║\n');
    fprintf('╚══════════════════════════════════════════════════════════════╝\n\n');
    
    %% 加载参数扫描数据
    dataFile = fullfile(reportDir, 'cma_parameter_study_data.mat');
    if ~exist(dataFile, 'file')
        error('请先运行 generate_cma_parameter_study.m 生成数据');
    end
    
    S = load(dataFile);
    data = S.report_data;
    
    %% ==================== 创建简化分析报告图 ====================
    
    figure('Name', 'CMA参数深度分析', 'Position', [50, 50, 1400, 800]);
    
    %% ========== 1. 步长分析 ==========
    step_study = data.step_study;
    step_values = step_study.values;
    step_labels = step_study.labels;
    
    % 过滤掉NaN值用于分析
    valid_step = ~isnan(step_study.evm_after);
    
    % 1.1 步长-EVM-收敛
    subplot(2, 3, 1);
    yyaxis left;
    semilogx(step_values(valid_step), step_study.evm_after(valid_step), 'b-o', 'LineWidth', 2, 'MarkerSize', 8);
    ylabel('EVM (%)', 'Color', 'b');
    ylim([45 47]);
    
    yyaxis right;
    semilogx(step_values(valid_step), step_study.steady_state_err(valid_step), 'r-s', 'LineWidth', 2, 'MarkerSize', 8);
    ylabel('稳态误差', 'Color', 'r');
    
    xlabel('步长 μ');
    title('步长 vs EVM/稳态误差');
    grid on;
    
    % 1.2 不同步长收敛曲线
    subplot(2, 3, 2);
    colors = lines(sum(valid_step));
    color_idx = 1;
    legend_entries = {};
    for i = find(valid_step)'
        if ismember(i, [1, 3, 5, 6])
            result = step_study.results{i};
            err_smooth = movmean(abs(result.errx_vec), 2000);
            plot((1:length(err_smooth))/1000, err_smooth, 'LineWidth', 1.5, 'Color', colors(color_idx,:));
            hold on;
            legend_entries{end+1} = step_labels{i};
            color_idx = color_idx + 1;
        end
    end
    xlabel('样本索引 (×10³)');
    ylabel('|误差|');
    title('不同步长收敛曲线');
    legend(legend_entries, 'Location', 'northeast');
    grid on;
    xlim([0 100]);
    
    %% ========== 2. 抽头长度分析 ==========
    tap_study = data.tap_study;
    tap_values = tap_study.values;
    valid_tap = ~isnan(tap_study.evm_after);
    
    % 2.1 抽头长度-EVM-时间
    subplot(2, 3, 3);
    yyaxis left;
    plot(tap_values(valid_tap), tap_study.evm_after(valid_tap), 'b-o', 'LineWidth', 2, 'MarkerSize', 8);
    ylabel('EVM (%)', 'Color', 'b');
    
    yyaxis right;
    plot(tap_values(valid_tap), tap_study.exec_time(valid_tap), 'r-s', 'LineWidth', 2, 'MarkerSize', 8);
    ylabel('执行时间 (s)', 'Color', 'r');
    
    xlabel('抽头长度 L');
    title('抽头长度 vs EVM/时间');
    grid on;
    
    %% ========== 3. 分段长度与并行段数分析 ==========
    seg_study = data.seg_study;
    par_study = data.par_study;
    
    seg_values = seg_study.values;
    par_values = par_study.values;
    valid_seg = ~isnan(seg_study.evm_after);
    valid_par = ~isnan(par_study.evm_after);
    
    % 3.1 分段长度分析
    subplot(2, 3, 4);
    yyaxis left;
    semilogx(seg_values(valid_seg), seg_study.evm_after(valid_seg), 'b-o', 'LineWidth', 2, 'MarkerSize', 8);
    ylabel('EVM (%)', 'Color', 'b');
    ylim([45 46]);
    
    yyaxis right;
    semilogx(seg_values(valid_seg), seg_study.exec_time(valid_seg), 'r-s', 'LineWidth', 2, 'MarkerSize', 8);
    ylabel('执行时间 (s)', 'Color', 'r');
    
    xlabel('分段长度 B');
    title('分段长度 vs EVM/时间');
    grid on;
    set(gca, 'XTick', seg_values);
    
    % 3.2 并行段数分析
    subplot(2, 3, 5);
    yyaxis left;
    semilogx(par_values(valid_par), par_study.evm_after(valid_par), 'b-o', 'LineWidth', 2, 'MarkerSize', 8);
    ylabel('EVM (%)', 'Color', 'b');
    
    yyaxis right;
    semilogx(par_values(valid_par), par_study.convergence_idx(valid_par)/1000, 'r-s', 'LineWidth', 2, 'MarkerSize', 8);
    ylabel('收敛点 (×10³)', 'Color', 'r');
    
    xlabel('并行段数 P');
    title('并行段数 vs EVM/收敛');
    grid on;
    set(gca, 'XTick', par_values);
    
    % 3.3 抽头更新周期分析
    subplot(2, 3, 6);
    update_period = 32 * par_values(valid_par);
    plot(update_period, par_study.convergence_idx(valid_par)/1000, 'g-^', 'LineWidth', 2, 'MarkerSize', 8);
    xlabel('抽头更新周期 T=B×P (样本)');
    ylabel('收敛点 (×10³)');
    title('更新周期 vs 收敛速度');
    grid on;
    
    sgtitle('CMA均衡器参数性能分析', 'FontSize', 14, 'FontWeight', 'bold');
    
    savefig(fullfile(reportDir, 'cma_parameter_deep_analysis.fig'));
    print(fullfile(reportDir, 'cma_parameter_deep_analysis.png'), '-dpng', '-r150');
    fprintf('深度分析图保存到: %s\n', fullfile(reportDir, 'cma_parameter_deep_analysis.png'));
    
    %% ==================== 打印详细分析报告 ====================
    fprintf('\n');
    fprintf('╔══════════════════════════════════════════════════════════════════════════════╗\n');
    fprintf('║                     CMA均衡器参数性能深度分析报告                            ║\n');
    fprintf('╚══════════════════════════════════════════════════════════════════════════════╝\n\n');
    
    %% 步长分析
    fprintf('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');
    fprintf('【1. 步长 (μ) 性能分析】\n');
    fprintf('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n');
    
    fprintf('■ 算法原理:\n');
    fprintf('  CMA抽头更新公式: w(n+1) = w(n) + μ · e(n) · y(n) · x*(n)\n');
    fprintf('  其中 e(n) = R - |y(n)|² 是恒模误差\n\n');
    
    fprintf('■ 测试结果:\n');
    fprintf('  ┌─────────┬──────────┬────────────┬────────────┐\n');
    fprintf('  │  步长   │  EVM(%%)  │ 收敛点(k)  │ 稳态误差   │\n');
    fprintf('  ├─────────┼──────────┼────────────┼────────────┤\n');
    for i = 1:length(step_values)
        if valid_step(i)
            fprintf('  │ %-7s │  %5.2f   │   %6.1f   │   %.4f   │\n', ...
                    step_labels{i}, step_study.evm_after(i), ...
                    step_study.convergence_idx(i)/1000, step_study.steady_state_err(i));
        else
            fprintf('  │ %-7s │   NaN    │     发散   │     -      │\n', step_labels{i});
        end
    end
    fprintf('  └─────────┴──────────┴────────────┴────────────┘\n\n');
    
    fprintf('■ 趋势成因分析:\n\n');
    fprintf('  1) 步长与收敛速度的关系:\n');
    fprintf('     理论: 收敛时间常数 τ ≈ 1/(μ·λmin), λmin为输入协方差矩阵最小特征值\n');
    fprintf('     实测: μ从2⁻¹²增至2⁻⁸, 收敛点从33k降至19k (加速1.7倍)\n');
    fprintf('     原因: 步长增大使抽头调整幅度增大, 更快接近最优解\n\n');
    
    fprintf('  2) 步长与稳态误差的关系:\n');
    fprintf('     理论: 稳态MSE ≈ μ·σ²n·trace(R)/4, 与步长成正比\n');
    fprintf('     实测: μ从2⁻¹²增至2⁻⁷, 稳态误差从0.308增至0.389 (增加26%%)\n');
    fprintf('     原因: 大步长导致抽头在最优点附近振荡, 无法精确收敛\n\n');
    
    fprintf('  3) 步长过大导致发散:\n');
    fprintf('     理论: 稳定条件 μ < 2/(λmax), λmax为最大特征值\n');
    fprintf('     实测: μ≥2⁻⁶时算法发散(输出NaN)\n');
    fprintf('     原因: 抽头更新过度补偿, 误差不减反增, 形成正反馈导致发散\n\n');
    
    fprintf('■ 最优选择: μ = 2⁻¹⁰ ~ 2⁻⁸\n');
    fprintf('  在此区间内, EVM≈45.35%%, 收敛时间<20k样本, 实现速度-精度平衡\n\n');
    
    %% 抽头长度分析
    fprintf('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');
    fprintf('【2. 抽头长度 (L) 性能分析】\n');
    fprintf('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n');
    
    fprintf('■ 物理意义:\n');
    fprintf('  均衡器时域跨度 = L × Ts (Ts为符号周期)\n');
    fprintf('  对于28GBaud系统: L=33 → 覆盖约1.2ns时延扩展\n\n');
    
    fprintf('■ 测试结果:\n');
    fprintf('  ┌──────────┬──────────┬────────────┬────────────┐\n');
    fprintf('  │ 抽头长度 │  EVM(%%)  │ 执行时间(s)│  收敛点(k) │\n');
    fprintf('  ├──────────┼──────────┼────────────┼────────────┤\n');
    for i = 1:length(tap_values)
        if valid_tap(i)
            fprintf('  │    %3d   │  %5.2f   │   %6.3f   │   %6.1f   │\n', ...
                    tap_values(i), tap_study.evm_after(i), ...
                    tap_study.exec_time(i), tap_study.convergence_idx(i)/1000);
        else
            fprintf('  │    %3d   │   NaN    │   %6.3f   │    发散    │\n', ...
                    tap_values(i), tap_study.exec_time(i));
        end
    end
    fprintf('  └──────────┴──────────┴────────────┴────────────┘\n\n');
    
    fprintf('■ 趋势成因分析:\n\n');
    fprintf('  1) 抽头过短(L=9)时EVM略高:\n');
    fprintf('     原因: 均衡器覆盖范围不足, 无法完全补偿信道的码间干扰(ISI)\n');
    fprintf('     表现: EVM=45.46%%, 比最优高0.08%%\n\n');
    
    fprintf('  2) 抽头在17-33范围最优:\n');
    fprintf('     原因: 该长度足以覆盖光纤信道的有效时延扩展\n');
    fprintf('     表现: EVM稳定在45.38-45.41%%\n\n');
    
    fprintf('  3) 抽头过长(≥49)性能下降:\n');
    fprintf('     原因: 参数过多导致:\n');
    fprintf('           a) 过拟合噪声而非信道\n');
    fprintf('           b) 收敛所需样本数增加\n');
    fprintf('           c) 梯度估计方差增大\n');
    fprintf('     表现: L=97时EVM升至46.07%%, L=129时发散\n\n');
    
    fprintf('  4) 计算时间非单调:\n');
    fprintf('     L=9时间最长(34s): FFT长度未优化, 非2的幂次\n');
    fprintf('     L=33时间最短(3.6s): FFT长度为64, 效率最高\n');
    fprintf('     L增大时间增加: O(L·log L) 复杂度\n\n');
    
    fprintf('■ 最优选择: L = 17~33\n');
    fprintf('  建议选择33: 最佳EVM且计算最快(FFT长度=64为2的幂)\n\n');
    
    %% 分段长度与并行段数分析
    fprintf('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');
    fprintf('【3. 分段长度 (B) 与 并行段数 (P) 分析】\n');
    fprintf('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n');
    
    fprintf('■ 关键关系: 抽头更新周期 T = B × P (样本)\n\n');
    
    fprintf('■ 分段长度测试结果:\n');
    fprintf('  ┌──────────┬──────────┬────────────┬──────────────┐\n');
    fprintf('  │ 分段长度 │  EVM(%%)  │ 执行时间(s)│ 更新周期(B×4)│\n');
    fprintf('  ├──────────┼──────────┼────────────┼──────────────┤\n');
    for i = 1:length(seg_values)
        if valid_seg(i)
            fprintf('  │    %3d   │  %5.2f   │   %6.3f   │     %4d     │\n', ...
                    seg_values(i), seg_study.evm_after(i), ...
                    seg_study.exec_time(i), seg_values(i)*4);
        else
            fprintf('  │    %3d   │   NaN    │   %6.3f   │     %4d     │\n', ...
                    seg_values(i), seg_study.exec_time(i), seg_values(i)*4);
        end
    end
    fprintf('  └──────────┴──────────┴────────────┴──────────────┘\n\n');
    
    fprintf('■ 并行段数测试结果:\n');
    fprintf('  ┌──────────┬──────────┬────────────┬───────────────┐\n');
    fprintf('  │ 并行段数 │  EVM(%%)  │ 收敛点(k)  │ 更新周期(32×P)│\n');
    fprintf('  ├──────────┼──────────┼────────────┼───────────────┤\n');
    for i = 1:length(par_values)
        if valid_par(i)
            fprintf('  │    %2d    │  %5.2f   │   %6.1f   │      %4d     │\n', ...
                    par_values(i), par_study.evm_after(i), ...
                    par_study.convergence_idx(i)/1000, 32*par_values(i));
        else
            fprintf('  │    %2d    │   NaN    │    发散    │      %4d     │\n', ...
                    par_values(i), 32*par_values(i));
        end
    end
    fprintf('  └──────────┴──────────┴────────────┴───────────────┘\n\n');
    
    fprintf('■ 趋势成因分析:\n\n');
    fprintf('  1) 分段长度B影响:\n');
    fprintf('     B太小(8): FFT开销比例高, 执行慢(10.5s)\n');
    fprintf('     B适中(32-64): 效率最优, 时间2-4s\n');
    fprintf('     B太大(≥128): 更新周期T≥512, 无法跟踪信道变化→发散\n\n');
    
    fprintf('  2) 并行段数P影响:\n');
    fprintf('     P=1时发散: 单段梯度噪声过大, 更新不稳定\n');
    fprintf('     P=2~8: 梯度平均降低噪声, 收敛稳定\n');
    fprintf('     P增大: 收敛变慢(P=32时收敛点24k vs P=4时19k)\n\n');
    
    fprintf('  3) 发散机理:\n');
    fprintf('     抽头更新周期T过大时, 信道可能在两次更新间发生显著变化\n');
    fprintf('     均衡器无法跟踪, 误差累积导致发散\n');
    fprintf('     临界周期: T ≈ 512 (B=128,P=4 或 B=32,P=16 边界)\n\n');
    
    fprintf('■ 最优选择: B=32, P=4~8\n');
    fprintf('  更新周期T=128~256, 在稳定性与跟踪能力间取得平衡\n\n');
    
    %% 综合建议
    fprintf('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');
    fprintf('【4. 综合参数调优建议】\n');
    fprintf('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n');
    
    fprintf('  ┌────────────┬─────────────┬────────────────────────────────────┐\n');
    fprintf('  │   参数     │  推荐值     │  调优原则                          │\n');
    fprintf('  ├────────────┼─────────────┼────────────────────────────────────┤\n');
    fprintf('  │  步长 μ    │  2⁻⁹~2⁻⁸    │  先大后小: 快速收敛后减小步长      │\n');
    fprintf('  │  抽头 L    │   33        │  ≥信道长度, 选2ⁿ-1以优化FFT        │\n');
    fprintf('  │  分段 B    │   32~64     │  权衡延迟与效率, 选2ⁿ              │\n');
    fprintf('  │  并行 P    │   4~8       │  P≥2确保稳定, P≤8保证跟踪能力      │\n');
    fprintf('  └────────────┴─────────────┴────────────────────────────────────┘\n\n');
    
    fprintf('  ★ 快速收敛场景: μ=2⁻⁸, L=33, B=32, P=4\n');
    fprintf('  ★ 高精度场景:   μ=2⁻¹⁰, L=33, B=32, P=8\n');
    fprintf('  ★ 低延迟场景:   μ=2⁻⁸, L=17, B=16, P=2\n\n');
    
    fprintf('═════════════════════════════════════════════════════════════════════════════\n');
    fprintf('                              分析报告完成                                   \n');
    fprintf('═════════════════════════════════════════════════════════════════════════════\n\n');

end

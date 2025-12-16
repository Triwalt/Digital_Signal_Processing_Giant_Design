% test_fft.m: FFT/IFFT 正确性与性能综合验证脚本
%
% 【目标】
%   除了“是否通过”，更关注“误差统计特性 + 可视化重合展示 + 性能/缓存影响”。
%   本脚本包含两部分：
%   1) 基2 FFT (my_fft) 与 MATLAB fft 的一致性验证
%   2) 任意点数 FFT (my_fft_mix, legacy) 与 MATLAB fft 的一致性验证
%
% 【验证维度】
%   A. 误差统计: max / RMS / mean / std / 95%分位数、相对误差等
%   B. 谱一致性: 幅度(dB)与相位的抽样点重合对比 + 误差曲线
%   C. 性能: 计时对比（包含 warmup 与缓存影响）
%   D. 一致性恒等式: Parseval 能量守恒、IFFT 往返重构
%
% 【说明】
%   - 性能测试受 MATLAB 版本、JIT、CPU、负载影响较大；主要用于相对对比。
%   - 绘图仅对部分长度/信号进行抽样展示，避免图太多。

clear; clc;
fprintf('========================================\n');
fprintf('FFT/IFFT 正确性与性能综合验证\n');
fprintf('========================================\n\n');

% =====================================================================
% 测试配置
% =====================================================================
% 基2 FFT 测试点数 (必须 2^m)
test_lengths_pow2 = [64, 256, 1024, 4096];

% 任意点数 FFT 测试点数 (包含题主关心的 N=80)
test_lengths_mixed = [15, 80, 96, 100, 255, 257];

% 信号种类（用于统计）
signal_kinds = { ...
    'impulse', ...
    'dc', ...
    'single_tone', ...
    'two_tone', ...
    'chirp_like', ...
    'white_noise' ...
};

% 误差阈值（仅作为“是否异常”的参考，统计更重要）
% 对于大 N 值（如 4096），浮点误差累积可能达到 1e-12 级别
tolerance_abs = 1e-10;

% 性能测试配置
bench_repeats = 30;   % 计时重复次数
bench_warmup = 3;     % 预热次数（触发 JIT、缓存构建等）

% 绘图配置（只对少数 case 展示，避免生成大量 figure）
do_plot = true;
plot_cases = [ ...
    struct('N', 1024, 'kind', 'white_noise', 'algo', 'my_fft'); ...
    struct('N', 80,   'kind', 'white_noise', 'algo', 'my_fft_mix') ...
];

all_passed = true;

% 汇总结果容器
results = struct('algo', {}, 'N', {}, 'kind', {}, 'stats', {}, 'timing', {});

% =====================================================================
% 工具函数 (脚本内局部函数)
% =====================================================================
make_signal = @local_make_signal;
err_stats = @local_error_stats;
bench = @local_benchmark;
report_stats = @local_report_stats;
plot_compare = @local_plot_compare;

% =====================================================================
% Part 1: 基2 FFT (my_fft)
% =====================================================================
fprintf('【Part 1】基2 FFT (my_fft) vs MATLAB fft\n');
fprintf('--------------------------------------------------------\n');

for N = test_lengths_pow2
    fprintf('\n--- N = %d (2^m) ---\n', N);

    for si = 1:numel(signal_kinds)
        kind = signal_kinds{si};
        x = make_signal(N, kind);

        Y_custom = my_fft(x);
        Y_builtin = fft(x);

        s = err_stats(Y_custom, Y_builtin);
        report_stats('my_fft', N, kind, s);

        % 误差门限：仅作为“异常告警”
        if s.max_abs > tolerance_abs
            all_passed = false;
        end

        % Parseval 能量守恒: sum|x|^2 ≈ (1/N)sum|X|^2
        ex = sum(abs(x).^2);
        eX = sum(abs(Y_custom).^2) / N;
        parseval_err = abs(eX - ex) / max(ex, eps);
        if parseval_err > 1e-10
            fprintf('    Parseval 警告: rel_err=%.2e\n', parseval_err);
        end

        % IFFT 往返（仅对部分类型做，避免过慢）
        if any(strcmp(kind, {'white_noise', 'chirp_like'}))
            x_rec = my_ifft(Y_custom);
            rt_err = max(abs(x_rec - x));
            if rt_err > 1e-10
                fprintf('    IFFT 往返警告: max_abs=%.2e\n', rt_err);
            end
        end

        % 记录结果
        results(end+1) = struct('algo','my_fft','N',N,'kind',kind,'stats',s,'timing',[]); %#ok<SAGROW>
    end

    % 性能基准（基2 FFT）
    t_my = bench(@() my_fft(make_signal(N, 'white_noise')), bench_warmup, bench_repeats);
    t_fft = bench(@() fft(make_signal(N, 'white_noise')), bench_warmup, bench_repeats);
    fprintf('  [Timing] my\\_fft median=%.3g ms | fft median=%.3g ms | ratio=%.2f\n', ...
        1e3*t_my.median, 1e3*t_fft.median, t_my.median / t_fft.median);

    % 缓存效应（my_fft 有 persistent cache，第一轮可能更慢）
    t_first = timeit(@() my_fft(make_signal(N, 'white_noise')));
    t_second = timeit(@() my_fft(make_signal(N, 'white_noise')));
    fprintf('  [Cache ] my\\_fft timeit first=%.3g ms | second=%.3g ms\n', 1e3*t_first, 1e3*t_second);
end

% =====================================================================
% Part 2: 任意点数 FFT (my_fft_mix)
% =====================================================================
fprintf('\n【Part 2】任意点数 FFT (my_fft_mix, legacy) vs MATLAB fft\n');
fprintf('--------------------------------------------------------\n');

have_mix = exist('my_fft_mix', 'file') == 2;
if ~have_mix
    fprintf('警告: 未找到 my_fft_mix.m，跳过 Part 2\n');
else
    for N = test_lengths_mixed
        fprintf('\n--- N = %d (mixed) ---\n', N);

        for si = 1:numel(signal_kinds)
            kind = signal_kinds{si};
            x = make_signal(N, kind);

            Y_custom = my_fft_mix(x);
            Y_builtin = fft(x);

            s = err_stats(Y_custom, Y_builtin);
            report_stats('my_fft_mix', N, kind, s);

            % Bluestein 算法（用于大质数）会有更多浮点运算，
            % 误差可能达到 1e-11 级别
            if s.max_abs > 1e-8
                all_passed = false;
            end

            results(end+1) = struct('algo','my_fft_mix','N',N,'kind',kind,'stats',s,'timing',[]); %#ok<SAGROW>
        end

        % 性能基准（注意 my_fft_mix 递归 + exp + 可能 Bluestein，常数较大）
        t_my = bench(@() my_fft_mix(make_signal(N, 'white_noise')), bench_warmup, bench_repeats);
        t_fft = bench(@() fft(make_signal(N, 'white_noise')), bench_warmup, bench_repeats);
        fprintf('  [Timing] my\\_fft\\_mix median=%.3g ms | fft median=%.3g ms | ratio=%.2f\n', ...
            1e3*t_my.median, 1e3*t_fft.median, t_my.median / t_fft.median);
    end
end

% =====================================================================
% 可视化对比
% =====================================================================
if do_plot
    for ci = 1:numel(plot_cases)
        pc = plot_cases(ci);
        N = pc.N;
        kind = pc.kind;
        x = make_signal(N, kind);
        Y_builtin = fft(x);
        if strcmp(pc.algo, 'my_fft')
            if mod(log2(N), 1) ~= 0
                continue;
            end
            Y_custom = my_fft(x);
        else
            if exist('my_fft_mix', 'file') ~= 2
                continue;
            end
            Y_custom = my_fft_mix(x);
        end
        plot_compare(pc.algo, N, kind, Y_custom, Y_builtin);
    end
end

% 总结
fprintf('========================================\n');
if all_passed
    fprintf('所有测试未发现明显异常 ✓\n');
else
    fprintf('存在超阈值误差/警告项 ✗ (建议查看上方统计输出)\n');
end
fprintf('========================================\n');


% =====================================================================
% Local functions
% =====================================================================
function x = local_make_signal(N, kind)
% 统一输出列向量 (与 my_fft 风格一致)
rng(42);
n = (0:N-1).';

switch kind
    case 'impulse'
        x = zeros(N, 1);
        x(1) = 1;
    case 'dc'
        x = ones(N, 1);
    case 'single_tone'
        k = max(1, floor(N/8));
        x = exp(1j * 2 * pi * k * n / N);
    case 'two_tone'
        k1 = max(1, floor(N/10));
        k2 = max(1, floor(N/6));
        x = exp(1j * 2 * pi * k1 * n / N) + 0.7 * exp(1j * 2 * pi * k2 * n / N);
    case 'chirp_like'
        % 简单“啁啾”相位（不依赖 Signal Processing Toolbox）
        x = exp(1j * 2 * pi * (0.2 * (n.^2) / max(N,1)));
    case 'white_noise'
        x = randn(N, 1) + 1j * randn(N, 1);
    otherwise
        error('Unknown signal kind: %s', kind);
end
end


function s = local_error_stats(Y_custom, Y_ref)
e = Y_custom - Y_ref;
abs_e = abs(e);
abs_ref = abs(Y_ref);

s = struct();
s.max_abs = max(abs_e);
s.rms = sqrt(mean(abs_e.^2));
s.mean_abs = mean(abs_e);
s.std_abs = std(abs_e);
s.p95_abs = prctile(abs_e, 95);

% 相对误差（对 ref 很小的点做下限保护）
den = max(abs_ref, 1e-15);
rel = abs_e ./ den;
s.max_rel = max(rel);
s.rms_rel = sqrt(mean(rel.^2));

% 复误差的均值（看是否存在系统性偏差）
s.mean_complex = mean(e);
end


function local_report_stats(algo, N, kind, s)
fprintf('  %-10s | %-11s | max=%.2e  rms=%.2e  mean=%.2e  std=%.2e  p95=%.2e  maxRel=%.2e\n', ...
    algo, kind, s.max_abs, s.rms, s.mean_abs, s.std_abs, s.p95_abs, s.max_rel);
end


function t = local_benchmark(fn, warmup, repeats)
for i = 1:warmup
    fn();
end

times = zeros(repeats, 1);
for i = 1:repeats
    tic;
    fn();
    times(i) = toc;
end

t = struct();
t.times = times;
t.median = median(times);
t.mean = mean(times);
t.std = std(times);
end


function local_plot_compare(algo, N, kind, Y_custom, Y_ref)
e = Y_custom - Y_ref;

% 抽样点：对数间隔 + 端点，保证小 N/大 N 都有代表性
numPts = min(200, N);
idx = unique(round(logspace(0, log10(N), numPts)));
idx(idx < 1) = [];
idx(idx > N) = [];

mag_ref = abs(Y_ref);
mag_cus = abs(Y_custom);
mag_db_ref = 20*log10(mag_ref + 1e-15);
mag_db_cus = 20*log10(mag_cus + 1e-15);
phase_ref = angle(Y_ref);
phase_cus = angle(Y_custom);

figure('Name', sprintf('%s N=%d %s', algo, N, kind), 'Color', 'w');
tiledlayout(3,1,'TileSpacing','compact','Padding','compact');

nexttile;
plot(idx, mag_db_ref(idx), 'k-', 'LineWidth', 1.0); hold on;
plot(idx, mag_db_cus(idx), 'r--', 'LineWidth', 1.0);
grid on; ylabel('|X| (dB)');
title(sprintf('%s vs fft | N=%d | %s (sampled points)', algo, N, kind));
legend('fft','custom','Location','best');

nexttile;
plot(idx, phase_ref(idx), 'k-', 'LineWidth', 1.0); hold on;
plot(idx, phase_cus(idx), 'r--', 'LineWidth', 1.0);
grid on; ylabel('phase (rad)');

nexttile;
semilogy(1:N, abs(e)+1e-18, 'b-');
grid on; xlabel('k'); ylabel('|err|');
end

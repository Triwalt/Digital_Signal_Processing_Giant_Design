%% cma_equalizer_core.m: 可配置参数的CMA均衡器核心函数
%
% 功能说明:
%   对双极化信号进行2x2 MIMO CMA盲均衡
%   支持自定义步长、抽头长度等关键参数
%
% 输入:
%   SigX_in - X极化输入信号 (列向量)
%   SigY_in - Y极化输入信号 (列向量)
%   params  - 参数结构体 (可选)
%       .step     - 步长 (默认 2^-8)
%       .tap_len  - 抽头长度 (默认 33, 必须为奇数)
%       .seg_len  - 分段长度 (默认 32)
%       .par_num  - 并行段数 (默认 4)
%       .clk_dly  - 时钟延迟 (默认 1)
%       .Rx       - X极化期望功率 (默认 2)
%       .Ry       - Y极化期望功率 (默认 2)
%       .use_builtin_fft - 使用MATLAB内置FFT (默认 true, 更快)
%
% 输出:
%   result - 结果结构体
%       .SigX_out   - X极化均衡输出
%       .SigY_out   - Y极化均衡输出
%       .errx_vec   - X极化误差向量
%       .erry_vec   - Y极化误差向量
%       .taps       - 最终抽头 (xx, yy, xy, yx)
%       .evm_before - 均衡前EVM
%       .evm_after  - 均衡后EVM
%       .convergence_idx - 收敛点索引
%       .exec_time  - 执行时间
%
% 作者: DSP课程设计项目
% 日期: 2025年12月

function result = cma_equalizer_core(SigX_in, SigY_in, params)

    %% 参数默认值
    if nargin < 3 || isempty(params)
        params = struct();
    end
    
    % 设置默认参数
    if ~isfield(params, 'step'), params.step = 2^-8; end
    if ~isfield(params, 'tap_len'), params.tap_len = 33; end
    if ~isfield(params, 'seg_len'), params.seg_len = 32; end
    if ~isfield(params, 'par_num'), params.par_num = 4; end
    if ~isfield(params, 'clk_dly'), params.clk_dly = 1; end
    if ~isfield(params, 'Rx'), params.Rx = 2; end
    if ~isfield(params, 'Ry'), params.Ry = 2; end
    if ~isfield(params, 'use_builtin_fft'), params.use_builtin_fft = true; end
    if ~isfield(params, 'verbose'), params.verbose = false; end
    
    % 参数提取
    step = params.step;
    tap_len = params.tap_len;
    seg_len = params.seg_len;
    par_num = params.par_num;
    clk_dly = params.clk_dly;
    Rx = params.Rx;
    Ry = params.Ry;
    
    % 确保tap_len为奇数
    if mod(tap_len, 2) == 0
        tap_len = tap_len + 1;
        if params.verbose
            warning('tap_len调整为奇数: %d', tap_len);
        end
    end
    
    tic;
    
    %% 数据预处理
    SigX_in = SigX_in(:);
    SigY_in = SigY_in(:);
    
    seg_num = floor(length(SigX_in) / seg_len);
    clk_num = ceil(seg_num / par_num);
    total_len = seg_len * seg_num;
    
    SigX_in = SigX_in(1:total_len);
    SigY_in = SigY_in(1:total_len);
    
    % 计算均衡前EVM
    constellation = [1+1j, 1-1j, -1-1j, -1+1j];  % 非归一化
    evm_before = calculate_evm_internal(SigX_in, constellation);
    
    % 前导零填充
    SigX_in_padded = [zeros(tap_len - 1, 1); SigX_in];
    SigY_in_padded = [zeros(tap_len - 1, 1); SigY_in];
    
    %% 初始化均衡器抽头
    xx = zeros(tap_len, clk_num + clk_dly);
    yy = zeros(tap_len, clk_num + clk_dly);
    xy = zeros(tap_len, clk_num + clk_dly);
    yx = zeros(tap_len, clk_num + clk_dly);
    
    center_tap = (tap_len + 1) / 2;
    xx(center_tap, :) = 1;
    yy(center_tap, :) = 1;
    
    %% 输出缓存
    SigX_out = zeros(total_len, 1);
    SigY_out = zeros(total_len, 1);
    errx_vec = zeros(total_len, 1);
    erry_vec = zeros(total_len, 1);
    
    %% 抽头增量缓存
    xx_accu = zeros(tap_len, 1);
    yy_accu = zeros(tap_len, 1);
    xy_accu = zeros(tap_len, 1);
    yx_accu = zeros(tap_len, 1);
    
    %% 频域卷积参数准备
    block_input_len = seg_len + tap_len - 1;
    conv_len = block_input_len + tap_len - 1;
    nfft_conv = 2^nextpow2(conv_len);
    
    % 零填充数组
    pad_input_x = zeros(nfft_conv, 1);
    pad_input_y = zeros(nfft_conv, 1);
    pad_xx = zeros(nfft_conv, 1);
    pad_xy = zeros(nfft_conv, 1);
    pad_yx = zeros(nfft_conv, 1);
    pad_yy = zeros(nfft_conv, 1);
    
    % 频域抽头缓存
    XX_fft_curr = [];
    XY_fft_curr = [];
    YX_fft_curr = [];
    YY_fft_curr = [];
    prev_clk_pos = -1;
    
    % 选择FFT函数
    if params.use_builtin_fft
        fft_func = @fft;
        ifft_func = @ifft;
    else
        fft_func = @my_fft;
        ifft_func = @my_ifft;
    end
    
    %% CMA主循环
    for k = 1:seg_num
        clk_pos = ceil(k / par_num);
        
        xx_curr = xx(:, clk_pos);
        yy_curr = yy(:, clk_pos);
        xy_curr = xy(:, clk_pos);
        yx_curr = yx(:, clk_pos);
        
        block_start = (k - 1) * seg_len + 1;
        block_end = block_start + seg_len - 1;
        
        % 提取输入数据块
        block_input_x = SigX_in_padded(block_start:block_end + tap_len - 1);
        block_input_y = SigY_in_padded(block_start:block_end + tap_len - 1);
        
        pad_input_x(:) = 0;
        pad_input_y(:) = 0;
        pad_input_x(1:block_input_len) = block_input_x;
        pad_input_y(1:block_input_len) = block_input_y;
        
        % 更新频域抽头缓存
        if clk_pos ~= prev_clk_pos
            pad_xx(:) = 0;
            pad_xy(:) = 0;
            pad_yx(:) = 0;
            pad_yy(:) = 0;
            
            pad_xx(1:tap_len) = xx_curr;
            pad_xy(1:tap_len) = xy_curr;
            pad_yx(1:tap_len) = yx_curr;
            pad_yy(1:tap_len) = yy_curr;
            
            XX_fft_curr = fft_func(pad_xx);
            XY_fft_curr = fft_func(pad_xy);
            YX_fft_curr = fft_func(pad_yx);
            YY_fft_curr = fft_func(pad_yy);
            prev_clk_pos = clk_pos;
        end
        
        % 频域卷积
        X_fft = fft_func(pad_input_x);
        Y_fft = fft_func(pad_input_y);
        
        conv_xx = ifft_func(X_fft .* XX_fft_curr);
        conv_xy = ifft_func(Y_fft .* XY_fft_curr);
        conv_yx = ifft_func(X_fft .* YX_fft_curr);
        conv_yy = ifft_func(Y_fft .* YY_fft_curr);
        
        % 截取有效输出
        block_out_x = conv_xx(tap_len:tap_len + seg_len - 1) + conv_xy(tap_len:tap_len + seg_len - 1);
        block_out_y = conv_yx(tap_len:tap_len + seg_len - 1) + conv_yy(tap_len:tap_len + seg_len - 1);
        
        SigX_out(block_start:block_end) = block_out_x;
        SigY_out(block_start:block_end) = block_out_y;
        
        % 计算误差
        errx_block = Rx - abs(block_out_x).^2;
        erry_block = Ry - abs(block_out_y).^2;
        
        errx_vec(block_start:block_end) = errx_block;
        erry_vec(block_start:block_end) = erry_block;
        
        % 计算梯度
        grad_x = step * (errx_block .* block_out_x);
        grad_y = step * (erry_block .* block_out_y);
        
        % 生成输入矩阵
        x_mat = toeplitz(block_input_x(tap_len:end), block_input_x(tap_len:-1:1));
        y_mat = toeplitz(block_input_y(tap_len:end), block_input_y(tap_len:-1:1));
        
        xx_accu = xx_accu + x_mat' * grad_x;
        xy_accu = xy_accu + y_mat' * grad_x;
        yx_accu = yx_accu + x_mat' * grad_y;
        yy_accu = yy_accu + y_mat' * grad_y;
        
        % 更新抽头
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
    end
    
    exec_time = toc;
    
    %% 计算均衡后EVM
    evm_after = calculate_evm_internal(SigX_out, constellation);
    
    %% 计算收敛点 (误差降到稳态的10%以内)
    err_smooth = movmean(abs(errx_vec), 1000);
    steady_state_err = mean(err_smooth(end-10000:end));
    convergence_threshold = steady_state_err * 1.1;
    convergence_idx = find(err_smooth < convergence_threshold, 1, 'first');
    if isempty(convergence_idx)
        convergence_idx = total_len;
    end
    
    %% 构建输出结构体
    result = struct();
    result.SigX_out = SigX_out;
    result.SigY_out = SigY_out;
    result.errx_vec = errx_vec;
    result.erry_vec = erry_vec;
    result.taps.xx = xx(:, clk_num);
    result.taps.yy = yy(:, clk_num);
    result.taps.xy = xy(:, clk_num);
    result.taps.yx = yx(:, clk_num);
    result.evm_before = evm_before;
    result.evm_after = evm_after;
    result.convergence_idx = convergence_idx;
    result.exec_time = exec_time;
    result.total_len = total_len;
    result.params = params;
    
end

%% 内部EVM计算函数
function evm = calculate_evm_internal(signal, constellation)
    signal = signal(:);
    % 判决到最近星座点
    distances = abs(signal - constellation);
    [~, idx] = min(distances, [], 2);
    decided = constellation(idx).';
    
    % 计算EVM
    error = signal - decided;
    evm = sqrt(mean(abs(error).^2)) / sqrt(mean(abs(constellation).^2)) * 100;
end

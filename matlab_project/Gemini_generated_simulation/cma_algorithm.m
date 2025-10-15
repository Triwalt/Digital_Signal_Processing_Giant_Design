function analysis = cma_algorithm(rx_signal, params)
% cma_algorithm: Implements the Constant Modulus Algorithm for blind equalization
%
% ALGORITHM PRINCIPLE:
% The Constant Modulus Algorithm (CMA) is a blind adaptive equalization technique
% that exploits the constant envelope property of certain modulated signals.
% It minimizes the cost function: J = E[(|y(n)|² - R²)²]
% where y(n) is the equalizer output and R² is the desired constant modulus.
%
% IMPLEMENTATION DETAILS:
% 1. CENTER-TAP INITIALIZATION: Equalizer weights initialized with center tap = 1
% 2. NORMALIZED CMA UPDATE: μ is normalized by input power per sample
% 3. PILOT/REFERENCE ASSISTED PHASE CORRECTION: least squares alignment
% 4. DETAILED DIAGNOSTICS: captures convergence, tap evolution, residual ISI
%
% INPUTS:
%   rx_signal    - Received signal vector (complex-valued, row or column)
%   params       - Structure with CMA configuration:
%                  .mu           - Base step size (μ)
%                  .num_taps     - Equalizer taps (>= channel length)
%                  .R2           - Constant modulus radius squared
%                  .guard        - Samples discarded before metric evaluation
%                  .max_iter     - Adaptation passes (default 1)
%                  .reference    - Optional known reference/pilot symbols
%                  .regularizer  - Small constant for numerical stability
%                  .weight_clip  - Maximum allowed tap magnitude
%
% OUTPUT (analysis struct):
%   .output        - Phase-aligned equalized signal after guard removal
%   .raw_output    - Raw CMA output before alignment
%   .weights       - Final equalizer weights
%   .error_signal  - CMA error sequence (post-guard)
%   .step_history  - Normalized step size evolution (post-guard)
%   .alpha         - Least-squares alignment coefficient
%   .evm           - RMS EVM (percentage)
%   .mse           - Residual CMA cost (mean |error|^2)
%   .decided       - Hard-decision constellation points
%   .guard         - Guard samples discarded
%   .params        - Echo of supplied parameters
%
% PERFORMANCE NOTES FROM LITERATURE:
% - Normalized CMA improves stability under power fluctuations
% - Weight clipping prevents tap explosion in fixed-point hardware
% - Pilot-aided phase rotation resolves CMA phase ambiguity for QAM

    arguments
        rx_signal (1,:) {mustBeNumeric}
        params.mu double {mustBePositive}
        params.num_taps (1,1) double {mustBeInteger,mustBePositive}
        params.R2 (1,1) double {mustBePositive} = 1
        params.guard (1,1) double {mustBeInteger} = 128
        params.max_iter (1,1) double {mustBeInteger,mustBePositive} = 1
        params.reference double = [] % optional reference symbols
        params.regularizer double {mustBePositive} = 1e-6
        params.weight_clip double {mustBePositive} = 10
    end

    % Ensure column vector
    if size(rx_signal, 1) == 1
        rx_signal = rx_signal(:);
    end
    N = length(rx_signal);
    L = params.num_taps;

    % Pad signal to handle initial conditions
    rx_padded = [zeros(L-1,1); rx_signal(:)];

    % Initialize weights (center spike)
    w = zeros(L,1);
    w(ceil(L/2)) = 1;

    % Preallocate
    y = zeros(N,1);
    err = zeros(N,1);
    step_history = zeros(N,1);

    % CMA adaptation (single pass default, extendable via max_iter)
    for iter = 1:params.max_iter
        for n = 1:N
            idx = n + L - 1;
            x_vec = rx_padded(idx:-1:idx-L+1);
            y_n = w' * x_vec;
            y(n) = y_n;

            % Normalized CMA step size (Eq. from literature): μ_n = μ / (||x||^2 + ε)
            power_norm = (x_vec' * x_vec) + params.regularizer;
            mu_n = params.mu / power_norm;

            % CMA error term (R2 is expected modulus)
            err_n = (params.R2 - abs(y_n)^2) * conj(y_n);
            err(n) = err_n;

            % Weight update using standard CMA梯度
            w = w + mu_n * err_n * x_vec;

            % Weight clipping to prevent divergence
            max_w = max(abs(w));
            if max_w > params.weight_clip
                w = w * (params.weight_clip / max_w);
            end

            step_history(n) = mu_n;
        end
    end

    % Discard guard transient from outputs
    guard = min(params.guard, floor(N/4));
    y_trunc = y(guard+1:end);
    err_trunc = err(guard+1:end);

    % Optional reference-based phase correction
    if ~isempty(params.reference)
        ref = params.reference(:);
        ref = ref(1:min(length(ref), length(y_trunc)));
        y_ref = y_trunc(1:length(ref));
        alpha = (y_ref' * ref) / (y_ref' * y_ref + params.regularizer);
    else
        alpha = exp(-1j*angle(mean(y_trunc(1:min(512,length(y_trunc))))));
    end
    y_aligned = alpha * y_trunc;

    % Constellation decision (QPSK assumption if reference absent)
    const_points = (1/sqrt(2)) * [1+1j, 1-1j, -1+1j, -1-1j];
    decided = zeros(size(y_aligned));
    for k = 1:length(y_aligned)
        [~, idx] = min(abs(y_aligned(k) - const_points.'));
        decided(k) = const_points(idx);
    end

    % Metrics
    evm_rms = sqrt(mean(abs(y_aligned - decided).^2)) / sqrt(mean(abs(decided).^2));
    mse = mean(abs(err_trunc).^2);

    analysis = struct();
    analysis.output = y_aligned;
    analysis.raw_output = y;
    analysis.truncated_output = y_trunc;
    analysis.weights = w;
    analysis.error_signal = err_trunc;
    analysis.step_history = step_history(guard+1:end);
    analysis.alpha = alpha;
    analysis.evm = evm_rms * 100;
    analysis.mse = mse;
    analysis.params = params;
    analysis.decided = decided;
end

function summary = generate_fft_vectors(varargin)
%GENERATE_FFT_VECTORS Produce fixed-point FFT stimuli and golden results.
%   SUMMARY = GENERATE_FFT_VECTORS('Name', Value, ...) creates stimulus and
%   reference data sets for parameterized FFT sizes. The script generates
%   multiple excitation scenarios, quantizes them to the requested fixed-
%   point format, computes the double-precision FFT as ground truth, then
%   quantizes the FFT outputs to match the hardware interface.
%
%   Optional name-value arguments:
%       'PointSizes'    : Vector of FFT lengths (default [32 128])
%       'WordLength'    : Total word length in bits (default 16)
%       'FracLength'    : Fractional length in bits (default 13 -> Q2.13)
%       'OutputDir'     : Directory for generated files
%                         (default '../data/fft_vectors')
%       'Scenarios'     : Cell array of custom scenario structs. When not
%                         provided, built-in scenarios are used.
%       'Overwrite'     : true/false to overwrite existing files (default true)
%
%   The function returns a SUMMARY struct containing metadata for all
%   generated vectors. Each scenario produces the following files per FFT
%   length N inside <OutputDir>/N/:
%       <scenario>_input_real.mem   -- Hexadecimal fixed-point real inputs
%       <scenario>_input_imag.mem
%       <scenario>_golden_real.mem  -- Hexadecimal fixed-point FFT outputs
%       <scenario>_golden_imag.mem
%       <scenario>_info.json        -- Metadata (scaling, RMS error, etc.)
%
%   Example:
%       generate_fft_vectors('PointSizes', [32 64], 'WordLength', 18, ...
%           'FracLength', 15, 'OutputDir', '../data/fft_vectors_q15');
%
%   See also FI, FFT.

% Copyright 2025

%% Parse arguments
p = inputParser;
p.addParameter('PointSizes', [32 128], @(x) isnumeric(x) && all(x == 2.^round(log2(x))));
p.addParameter('WordLength', 16, @(x) isnumeric(x) && isscalar(x) && x > 1);
p.addParameter('FracLength', 13, @(x) isnumeric(x) && isscalar(x) && x >= 0);
p.addParameter('OutputDir', fullfile(fileparts(mfilename('fullpath')), '..', 'data', 'fft_vectors'), @ischar);
p.addParameter('Scenarios', [], @(x) isempty(x) || iscell(x));
p.addParameter('Overwrite', true, @(x) islogical(x) && isscalar(x));
p.parse(varargin{:});
opts = p.Results;

if opts.FracLength >= opts.WordLength
    error('Fractional length must be smaller than total word length.');
end

%% Prepare directories
outputDir = opts.OutputDir;
if ~exist(outputDir, 'dir') %#ok<*EXIST>
    fprintf('[generate_fft_vectors] Creating output directory: %s\n', outputDir);
    mkdir(outputDir);
end

%% Fixed-point configuration
nt = numerictype('Signed', true, 'WordLength', opts.WordLength, ...
    'FractionLength', opts.FracLength);
fm = fimath('RoundingMethod', 'Nearest', 'OverflowAction', 'Saturate', ...
    'ProductMode', 'SpecifyPrecision', ...
    'ProductWordLength', opts.WordLength * 2, ...
    'ProductFractionLength', opts.FracLength * 2, ...
    'SumMode', 'SpecifyPrecision', ...
    'SumWordLength', opts.WordLength + 4, ...
    'SumFractionLength', opts.FracLength + 2);

%% Scenario configuration
if isempty(opts.Scenarios)
    scenarios = default_scenarios();
else
    scenarios = opts.Scenarios;
end

%% Summary container
summary = struct('wordLength', opts.WordLength, ...
                 'fracLength', opts.FracLength, ...
                 'pointSizes', sort(opts.PointSizes(:)'), ...
                 'scenarios', {scenarios}, ...
                 'outputs', []);

allOutputs = [];

%% Main generation loop
for nIdx = 1:numel(opts.PointSizes)
    N = opts.PointSizes(nIdx);
    nDir = fullfile(outputDir, sprintf('%d', N));
    if ~exist(nDir, 'dir')
        mkdir(nDir);
    end

    fprintf('[generate_fft_vectors] Processing N = %d\n', N);

    for sIdx = 1:numel(scenarios)
        scen = scenarios{sIdx};
        scenName = scen.name;
        fprintf('  -> Scenario: %s\n', scenName);

        % Generate complex stimulus (double precision)
        x = scen.generator(N);
        validateattributes(x, {'double'}, {'vector'});
        if numel(x) ~= N
            error('Scenario %s produced %d samples, expected %d.', scenName, numel(x), N);
        end
        x = x(:); % column vector

        % Fixed-point quantization
        [xinReal, xinImag, xq] = quantize_complex(x, nt, fm);

        % Golden FFT reference (double precision)
        Xgold = fft(double(xq));
        [xoutReal, xoutImag, Xq] = quantize_complex(Xgold, nt, fm);

        % Metrics
        err = Xgold - double(Xq);
        metrics = struct();
        metrics.maxAbsError = max(abs(err));
        metrics.rmsError = sqrt(mean(abs(err).^2));
        metrics.snr = 20 * log10(norm(double(Xgold)) / norm(err));
        metrics.scaling = struct('wordLength', opts.WordLength, ...
                                 'fracLength', opts.FracLength);

        % Paths
        baseName = sprintf('%s_%s', lower(scenName), scen.tag);
        realInPath = fullfile(nDir, sprintf('%s_input_real.mem', baseName));
        imagInPath = fullfile(nDir, sprintf('%s_input_imag.mem', baseName));
        realOutPath = fullfile(nDir, sprintf('%s_golden_real.mem', baseName));
        imagOutPath = fullfile(nDir, sprintf('%s_golden_imag.mem', baseName));
        infoPath = fullfile(nDir, sprintf('%s_info.json', baseName));

        if ~opts.Overwrite && exist(realInPath, 'file')
            warning('File %s exists. Skipping scenario.', realInPath);
            continue;
        end

        write_mem(realInPath, xinReal, opts.WordLength);
        write_mem(imagInPath, xinImag, opts.WordLength);
        write_mem(realOutPath, xoutReal, opts.WordLength);
        write_mem(imagOutPath, xoutImag, opts.WordLength);
        write_json(infoPath, scen, N, metrics, opts);

        allOutputs = [allOutputs; struct( ...
            'scenario', scenName, ...
            'tag', scen.tag, ...
            'N', N, ...
            'inputReal', realInPath, ...
            'inputImag', imagInPath, ...
            'goldenReal', realOutPath, ...
            'goldenImag', imagOutPath, ...
            'metrics', metrics)]; %#ok<AGROW>
    end
end

summary.outputs = allOutputs;

if nargout == 0
    clear summary;
end

fprintf('[generate_fft_vectors] Done. Generated %d dataset(s).\n', numel(allOutputs));

end

%% Local helpers ---------------------------------------------------------
function scenarios = default_scenarios()
scenarios = {
    struct('name', 'SingleTone', 'tag', 'tone1', 'generator', @single_tone), ...
    struct('name', 'DualTone',   'tag', 'tone2', 'generator', @dual_tone), ...
    struct('name', 'QAM4Random', 'tag', 'qam4',  'generator', @qam4_random), ...
    struct('name', 'Impulse',    'tag', 'imp',   'generator', @impulse), ...
    struct('name', 'Noise',      'tag', 'noise', 'generator', @white_noise) ...
    };
end

function [realOut, imagOut, fiComplex] = quantize_complex(x, nt, fm)
realFi = fi(real(x), nt, fm);
imagFi = fi(imag(x), nt, fm);
fiComplex = fi(realFi, nt, fm) + 1i * fi(imagFi, nt, fm);
realOut = storedInteger(realFi);
imagOut = storedInteger(imagFi);
end

function write_mem(filePath, data, wl)
fid = fopen(filePath, 'w');
if fid == -1
    error('Failed to open %s for writing.', filePath);
end
cleanup = onCleanup(@() fclose(fid));
width = wl / 4;
if mod(wl, 4) ~= 0
    error('Word length must be divisible by 4 for hex export.');
end
for k = 1:numel(data)
    word = uint32(typecast(int32(data(k)), 'uint32'));
    mask = bitshift(uint32(1), wl) - 1;
    word = bitand(word, mask);
    fprintf(fid, ['%0', num2str(width), 'X\n'], word);
end
fprintf('[generate_fft_vectors]   Wrote %s\n', filePath);
clear cleanup
end

function write_json(filePath, scen, N, metrics, opts)
info = struct();
info.timestamp = datestr(now, 'yyyy-mm-ddTHH:MM:SS');
info.pointSize = N;
info.scenario = scen.name;
info.tag = scen.tag;
info.generator = func2str(scen.generator);
info.wordLength = opts.WordLength;
info.fracLength = opts.FracLength;
info.metrics = metrics;
jsonText = jsonencode(info, 'PrettyPrint', true);
fid = fopen(filePath, 'w');
if fid == -1
    error('Failed to open %s for writing.', filePath);
end
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, '%s', jsonText);
fprintf('[generate_fft_vectors]   Wrote %s\n', filePath);
clear cleanup
end

%% Scenario generators ---------------------------------------------------
function x = single_tone(N)
% Pure tone at bin 4 (scaled to 0.8 full-scale)
n = (0:N-1).';
k = 4;
ampl = 0.8;
angles = 2 * pi * k * n / N;
x = ampl * exp(1i * angles);
end

function x = dual_tone(N)
% Sum of two tones at bins 3 and 9 with different amplitudes
n = (0:N-1).';
k1 = 3; k2 = 9;
ampl1 = 0.6;
ampl2 = 0.35;
x = ampl1 * exp(1i * 2 * pi * k1 * n / N) + ...
    ampl2 * exp(1i * 2 * pi * k2 * n / N + 1i * pi/4);
end

function x = qam4_random(N)
% Random 4-QAM symbols with Gray mapping, scaled to avoid clipping
M = 4;
symIdx = randi([0 M-1], N, 1);
map = (2 * floor(symIdx / 2) - 1) + 1i * (2 * mod(symIdx, 2) - 1);
x = 0.7 / sqrt(2) * map;
end

function x = impulse(N)
% Unit impulse followed by zeros
x = zeros(N, 1);
x(1) = 0.9;
end

function x = white_noise(N)
% Complex white Gaussian noise with 0.6 RMS
rng(2025); %#ok<RAND>
noise = (randn(N, 1) + 1i * randn(N, 1)) / sqrt(2);
x = 0.6 * noise;
end


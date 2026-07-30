function [rx_bits, ber, rx_data_symbols, debug] = ...
    wifi_rx_pipeline(rx_signal, cfg, tx_bits)
% WIFI_RX_PIPELINE
% OFDM receiver with two L-LTF preamble symbols followed by data symbols.
%
% Frame structure:
%   [L-LTF 1] [L-LTF 2] [DATA 1] ... [DATA cfg.numSymbols]
%
% The value cfg.numSymbols represents the number of DATA OFDM symbols.
% The two preamble symbols are added on top of cfg.numSymbols.
%
% Receiver tasks:
%   1. Validate the configuration and received frame length
%   2. Convert the complete serial frame to parallel OFDM symbols
%   3. Remove the cyclic prefix
%   4. Apply FFT and FFT shift
%   5. Reverse transmitter power normalisation
%   6. Separate the two L-LTF symbols from the data symbols
%   7. Estimate the channel from the two L-LTF symbols
%   8. Apply ZF or MMSE equalisation to DATA symbols only
%   9. Apply pilot-based residual phase correction
%  10. Extract data subcarriers
%  11. Perform QAM demodulation and bit recovery
%  12. Calculate BER
%
% Inputs:
%   rx_signal - Received serial time-domain frame
%   cfg       - Shared OFDM configuration structure
%   tx_bits   - Optional transmitted payload bits for BER calculation
%
% Outputs:
%   rx_bits         - Recovered payload bits
%   ber             - Bit error rate; NaN when tx_bits is not supplied
%   rx_data_symbols - Equalised data symbols before QAM demodulation
%   debug           - Intermediate receiver results
%
% Required cfg fields:
%   cfg.N_fft
%   cfg.N_cp
%   cfg.N_data
%   cfg.N_pilot
%   cfg.ModOrder
%   cfg.numSymbols       Number of DATA OFDM symbols
%   cfg.data_idx
%   cfg.pilot_idx
%
% Optional cfg fields used by equalise_ofdm.m:
%   cfg.equaliserType
%   cfg.noiseVarFreq
%   cfg.epsH
%   cfg.applyPhaseCorrection
%   cfg.pilot_base
%
% Required external functions:
%   estimate_channel.m
%   equalise_ofdm.m
%   generate_known_pilots.m
%
% Expected channel-estimator interface:
%   [Hhat, H_ltf] = estimate_channel(R_ltf, cfg)
%
% R_ltf has size:
%   cfg.N_fft x 2
%
% Hhat may have size:
%   cfg.N_fft x 1
% or
%   cfg.N_fft x cfg.numSymbols
%
% A single-column Hhat is automatically repeated across the data symbols.

    if nargin < 3
        tx_bits = [];
    end

    %% 1. Validate required configuration fields

    requiredFields = { ...
        'N_fft', ...
        'N_cp', ...
        'N_data', ...
        'N_pilot', ...
        'ModOrder', ...
        'numSymbols', ...
        'data_idx', ...
        'pilot_idx'};

    for fieldIndex = 1:length(requiredFields)

        fieldName = requiredFields{fieldIndex};

        if ~isfield(cfg, fieldName)
            error('Missing configuration field: cfg.%s', fieldName);
        end
    end

    if cfg.N_fft <= 0
        error('cfg.N_fft must be positive.');
    end

    if cfg.N_cp < 0
        error('cfg.N_cp must not be negative.');
    end

    if cfg.numSymbols <= 0 || mod(cfg.numSymbols, 1) ~= 0
        error('cfg.numSymbols must be a positive integer.');
    end

    if length(cfg.data_idx) ~= cfg.N_data
        error('Length of cfg.data_idx must equal cfg.N_data.');
    end

    if length(cfg.pilot_idx) ~= cfg.N_pilot
        error('Length of cfg.pilot_idx must equal cfg.N_pilot.');
    end

    if any(ismember(cfg.data_idx, cfg.pilot_idx))
        error('Data and pilot subcarrier indices must not overlap.');
    end

    if any(cfg.data_idx < 1) || any(cfg.data_idx > cfg.N_fft)
        error('cfg.data_idx contains an index outside 1:cfg.N_fft.');
    end

    if any(cfg.pilot_idx < 1) || any(cfg.pilot_idx > cfg.N_fft)
        error('cfg.pilot_idx contains an index outside 1:cfg.N_fft.');
    end

    bitsPerQamSymbol = log2(cfg.ModOrder);

    if mod(bitsPerQamSymbol, 1) ~= 0
        error('cfg.ModOrder must be a power of two.');
    end

    %% 2. Calculate complete frame dimensions
    %
    % cfg.numSymbols counts DATA symbols only.
    % The complete frame contains two additional L-LTF symbols.

    nPreamble = 2;
    totalSymbols = nPreamble + cfg.numSymbols;

    symbolLength = cfg.N_fft + cfg.N_cp;

    expectedLength = ...
        symbolLength * totalSymbols;

    %% 3. Check and truncate the received frame

    rx_signal = rx_signal(:);

    if length(rx_signal) < expectedLength

        error(['Received frame is too short. Expected at least %d ', ...
               'samples for %d preamble and %d data symbols, ', ...
               'but received %d samples.'], ...
               expectedLength, ...
               nPreamble, ...
               cfg.numSymbols, ...
               length(rx_signal));
    end

    % Ignore only samples beyond the complete 2-LTF-plus-data frame.
    rx_signal = rx_signal(1:expectedLength);

    %% 4. Serial-to-parallel conversion
    %
    % Matrix size:
    %   (N_fft + N_cp) x (2 + cfg.numSymbols)

    rx_frame_matrix = reshape( ...
        rx_signal, ...
        symbolLength, ...
        totalSymbols);

    %% 5. Cyclic-prefix removal

    rx_time_no_cp = ...
        rx_frame_matrix(cfg.N_cp+1:end, :);

    %% 6. FFT

    R_freq_unshifted = fft( ...
        rx_time_no_cp, ...
        cfg.N_fft, ...
        1);

    %% 7. FFT shift

    R_freq_shifted = fftshift( ...
        R_freq_unshifted, ...
        1);

    %% 8. Reverse transmitter power normalisation

    norm_factor = sqrt( ...
        cfg.N_fft / ...
        (cfg.N_data + cfg.N_pilot));

    R_freq_shifted = ...
        R_freq_shifted / norm_factor;

    R_freq_unshifted_scaled = ...
        R_freq_unshifted / norm_factor;

    %% 9. Separate L-LTF preamble and payload data symbols
    %
    % R_ltf:
    %   N_fft x 2
    %
    % R_data:
    %   N_fft x cfg.numSymbols

    R_ltf = ...
        R_freq_shifted(:, 1:nPreamble);

    R_data = ...
        R_freq_shifted(:, nPreamble+1:end);

    R_ltf_unshifted = ...
        R_freq_unshifted_scaled(:, 1:nPreamble);

    R_data_unshifted = ...
        R_freq_unshifted_scaled(:, nPreamble+1:end);

    if size(R_ltf, 2) ~= nPreamble
        error('The L-LTF matrix must contain exactly two columns.');
    end

    if size(R_data, 2) ~= cfg.numSymbols
        error(['The data grid contains %d symbols, but cfg.numSymbols ', ...
               'specifies %d symbols.'], ...
               size(R_data, 2), cfg.numSymbols);
    end

    %% 10. L-LTF-based channel estimation
    %
    % Only the two preamble columns are passed to the estimator.

    [Hhat_estimator, H_ltf] = ...
        estimate_channel( ...
            R_ltf, ...
            cfg);

    if size(Hhat_estimator, 1) ~= cfg.N_fft
        error('The channel estimate must have cfg.N_fft rows.');
    end

    %% 11. Match the channel-estimate size to the data grid

    if size(Hhat_estimator, 2) == 1

        % One packet-level channel estimate, reused for all data symbols.
        Hhat = repmat( ...
            Hhat_estimator, ...
            1, ...
            cfg.numSymbols);

    elseif size(Hhat_estimator, 2) == cfg.numSymbols

        % The estimator already returned one estimate per data symbol.
        Hhat = Hhat_estimator;

    else

        error(['The channel estimator returned %d columns. It must ', ...
               'return either one column or cfg.numSymbols (%d) columns.'], ...
               size(Hhat_estimator, 2), ...
               cfg.numSymbols);
    end

    if ~isequal(size(Hhat), size(R_data))
        error('The expanded channel estimate and data grid must have the same size.');
    end

    %% 12. Equalisation and residual phase correction
    %
    % Only payload data symbols are equalised.
    % The four scattered pilots are retained for per-symbol phase tracking.

    [R_final, X_zf_shifted, ...
     X_mmse_shifted, phase_error] = ...
        equalise_ofdm( ...
            R_data, ...
            Hhat, ...
            cfg);

    if size(R_final, 2) ~= cfg.numSymbols
        error('The equalised data grid has an incorrect number of symbols.');
    end

    %% 13. Extract the data subcarriers

    rx_data_grid = ...
        R_final(cfg.data_idx, :);

    %% 14. Parallel-to-serial conversion

    rx_data_symbols = ...
        rx_data_grid(:);

    %% 15. QAM demodulation and bit recovery

    rx_bits = qamdemod( ...
        rx_data_symbols, ...
        cfg.ModOrder, ...
        'OutputType', 'bit', ...
        'UnitAveragePower', true);

    rx_bits = double(rx_bits(:));

    %% 16. Preserve OFDM-symbol bit boundaries
    %
    % For QPSK:
    %   48 data subcarriers x 2 bits = 96 bits per OFDM symbol.

    bitsPerOFDMSymbol = ...
        cfg.N_data * bitsPerQamSymbol;

    expectedBitCount = ...
        bitsPerOFDMSymbol * cfg.numSymbols;

    if length(rx_bits) ~= expectedBitCount

        error(['Recovered %d bits, but %d payload bits were expected ', ...
               'from %d data OFDM symbols.'], ...
               length(rx_bits), ...
               expectedBitCount, ...
               cfg.numSymbols);
    end

    rx_bits_matrix = reshape( ...
        rx_bits, ...
        bitsPerOFDMSymbol, ...
        cfg.numSymbols);

    %% 17. BER calculation

    if isempty(tx_bits)

        ber = NaN;

    else

        tx_bits = double(tx_bits(:));

        if length(rx_bits) ~= length(tx_bits)

            error(['Recovered bit length (%d) does not match ', ...
                   'transmitted payload bit length (%d).'], ...
                   length(rx_bits), ...
                   length(tx_bits));
        end

        ber = mean(rx_bits ~= tx_bits);
    end

    %% 18. Debug outputs

    debug.nPreamble = nPreamble;
    debug.totalSymbols = totalSymbols;
    debug.symbolLength = symbolLength;
    debug.expectedLength = expectedLength;

    debug.rx_frame_matrix = rx_frame_matrix;
    debug.rx_time_no_cp = rx_time_no_cp;

    debug.R_freq_unshifted = ...
        R_freq_unshifted_scaled;

    debug.R_freq_shifted = ...
        R_freq_shifted;

    debug.R_ltf_unshifted = R_ltf_unshifted;
    debug.R_data_unshifted = R_data_unshifted;

    debug.R_ltf = R_ltf;
    debug.R_data = R_data;

    debug.Hhat_estimator = Hhat_estimator;
    debug.Hhat = Hhat;
    debug.H_ltf = H_ltf;

    debug.X_zf_shifted = X_zf_shifted;
    debug.X_mmse_shifted = X_mmse_shifted;

    debug.phase_error = phase_error;
    debug.R_final = R_final;

    debug.rx_data_grid = rx_data_grid;
    debug.rx_data_symbols = rx_data_symbols;

    debug.rx_bits_matrix = rx_bits_matrix;
    debug.bitsPerOFDMSymbol = bitsPerOFDMSymbol;

    debug.norm_factor = norm_factor;

end

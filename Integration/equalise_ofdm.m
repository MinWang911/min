function [R_final, X_zf_shifted, X_mmse_shifted, phase_error] = ...
    equalise_ofdm(R_freq_shifted, Hhat, cfg)
% EQUALISE_OFDM
% Apply ZF or MMSE equalisation to the received OFDM grid.
%
% Inputs:
%   R_freq_shifted : received OFDM frequency-domain grid after CP removal,
%                    FFT, fftshift, and power normalisation reversal
%   Hhat           : estimated channel on full fftshifted OFDM grid
%   cfg            : system configuration structure
%
% Required cfg fields:
%   cfg.N_fft
%   cfg.data_idx
%   cfg.pilot_idx
%
% Optional cfg fields:
%   cfg.equaliserType
%       'ZF' or 'MMSE'. Default is 'ZF'.
%
%   cfg.noiseVarFreq
%       Frequency-domain noise variance for MMSE.
%
%   cfg.epsH
%       Small value to avoid division by zero. Default is 1e-8.
%
%   cfg.applyPhaseCorrection
%       true or false. Default is true.
%
% Outputs:
%   R_final        : final equalised OFDM grid for receiver back-end
%   X_zf_shifted   : ZF equalised OFDM grid before phase correction
%   X_mmse_shifted : MMSE equalised OFDM grid before phase correction
%   phase_error    : residual phase error per OFDM symbol

    %% Part 1: Check required cfg fields

    requiredFields = {'N_fft', 'data_idx', 'pilot_idx'};

    for i = 1:length(requiredFields)
        if ~isfield(cfg, requiredFields{i})
            error(['Missing cfg field: ', requiredFields{i}]);
        end
    end

    Nfft = cfg.N_fft;
    data_idx = cfg.data_idx(:);
    pilot_idx = cfg.pilot_idx(:);

    if size(R_freq_shifted, 1) ~= Nfft
        error('The number of rows in R_freq_shifted must be equal to cfg.N_fft.');
    end

    if size(Hhat, 1) ~= Nfft
        error('The number of rows in Hhat must be equal to cfg.N_fft.');
    end

    if ~isequal(size(R_freq_shifted), size(Hhat))
        error('R_freq_shifted and Hhat must have the same size.');
    end

    numSymbols = size(R_freq_shifted, 2);

    %% Part 2: Read optional cfg fields

    if isfield(cfg, 'equaliserType')
        equaliserType = cfg.equaliserType;
    else
        equaliserType = 'ZF';
    end

    if isfield(cfg, 'noiseVarFreq')
        noiseVarFreq = cfg.noiseVarFreq;
    else
        noiseVarFreq = 0;
    end

    if isfield(cfg, 'epsH')
        epsH = cfg.epsH;
    else
        epsH = 1e-8;
    end

    if isfield(cfg, 'applyPhaseCorrection')
        applyPhaseCorrection = cfg.applyPhaseCorrection;
    else
        applyPhaseCorrection = true;
    end

    %% Part 3: Active subcarrier indices

    active_idx = sort([data_idx; pilot_idx]);

    %% Part 4: ZF equalisation

    H_safe = Hhat;

    H_safe(abs(H_safe) < epsH) = epsH;

    X_zf_shifted = zeros(Nfft, numSymbols);

    X_zf_shifted(active_idx, :) = ...
        R_freq_shifted(active_idx, :) ./ H_safe(active_idx, :);

    %% Part 5: MMSE equalisation

    H_active = Hhat(active_idx, :);
    Y_active = R_freq_shifted(active_idx, :);

    denominator = abs(H_active).^2 + noiseVarFreq;

    denominator(abs(denominator) < epsH) = epsH;

    X_mmse_shifted = zeros(Nfft, numSymbols);

    X_mmse_shifted(active_idx, :) = ...
        (conj(H_active) ./ denominator) .* Y_active;

    %% Part 6: Select ZF or MMSE output

    switch upper(equaliserType)

        case 'ZF'
            R_equalised = X_zf_shifted;

        case 'MMSE'
            R_equalised = X_mmse_shifted;

        otherwise
            error('cfg.equaliserType must be either ZF or MMSE.');

    end

    %% Part 7: Pilot-based residual phase correction

    if applyPhaseCorrection

        known_pilots = generate_known_pilots(cfg, numSymbols);

        equalised_pilots = R_equalised(pilot_idx, :);

        phase_error = angle( ...
            sum(equalised_pilots .* conj(known_pilots), 1));

        phase_correction = repmat(exp(-1j * phase_error), Nfft, 1);

        R_final = R_equalised .* phase_correction;

    else

        phase_error = zeros(1, numSymbols);

        R_final = R_equalised;

    end

end

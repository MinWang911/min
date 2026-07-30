function [Hhat, H_ltf] = estimate_channel(R_ltf_shifted, cfg)
% ESTIMATE_CHANNEL
% LTF-based channel estimation for OFDM receiver.
%
% This version replaces 4-pilot interpolation.
%
% Inputs:
%   R_ltf_shifted : received LTF frequency-domain grid after CP removal,
%                   FFT, fftshift, and normalisation reversal.
%                   Size: cfg.N_fft x 2
%
%   cfg           : system configuration structure
%
% Required cfg fields:
%   cfg.N_fft
%   cfg.data_idx
%   cfg.pilot_idx
%   cfg.numSymbols
%
% Outputs:
%   Hhat  : estimated channel on full fftshifted OFDM grid.
%           Size: cfg.N_fft x cfg.numSymbols
%
%   H_ltf : single channel estimate from the two LTF symbols.
%           Size: cfg.N_fft x 1

    %% Part 1: Check required fields

    requiredFields = {'N_fft', 'data_idx', 'pilot_idx', 'numSymbols'};

    for i = 1:length(requiredFields)
        if ~isfield(cfg, requiredFields{i})
            error(['Missing cfg field: ', requiredFields{i}]);
        end
    end

    Nfft = cfg.N_fft;

    data_idx = cfg.data_idx(:);
    pilot_idx = cfg.pilot_idx(:);

    active_idx = sort([data_idx; pilot_idx]);

    %% Part 2: Check LTF input

    if size(R_ltf_shifted, 1) ~= Nfft
        error('R_ltf_shifted must have cfg.N_fft rows.');
    end

    if size(R_ltf_shifted, 2) < 2
        error('R_ltf_shifted must contain two LTF symbols.');
    end

    %% Part 3: Generate known LTF sequence

    LTF_shifted = generate_known_ltf(cfg);

    %% Part 4: Extract received LTF symbols

    Y1 = R_ltf_shifted(:, 1);
    Y2 = R_ltf_shifted(:, 2);

    %% Part 5: LTF-based LS channel estimation

    H_ltf = zeros(Nfft, 1);

    H_ltf(active_idx) = ...
        (Y1(active_idx) + Y2(active_idx)) ./ ...
        (2 * LTF_shifted(active_idx));

    %% Part 6: Repeat H estimate for all data OFDM symbols

    Hhat = repmat(H_ltf, 1, cfg.numSymbols);

end

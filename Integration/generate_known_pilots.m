function known_pilots = generate_known_pilots(cfg, numSymbols)
% GENERATE_KNOWN_PILOTS
% Generate the known IEEE 802.11a pilot symbols.
%
% Inputs:
%   cfg        : system configuration structure
%   numSymbols : number of OFDM symbols
%
% Required cfg fields:
%   cfg.pilot_idx
%
% Optional cfg fields:
%   cfg.pilot_base
%
% Output:
%   known_pilots : known pilot matrix, size N_pilot x numSymbols

    if ~isfield(cfg, 'pilot_idx')
        error('Missing cfg field: pilot_idx');
    end

    if isfield(cfg, 'pilot_base')
        pilot_base = cfg.pilot_base(:);
    else
        pilot_base = [1; 1; 1; -1];
    end

    if length(pilot_base) ~= length(cfg.pilot_idx)
        error('Length of pilot_base must match length of cfg.pilot_idx.');
    end

    state = ones(1, 7);
    p_seq = zeros(1, numSymbols);

    for n = 1:numSymbols
        p_seq(n) = state(7);
        feedback = xor(state(7), state(4));
        state = [feedback, state(1:6)];
    end

    pilot_polarity = 1 - 2 * p_seq;

    known_pilots = pilot_base * pilot_polarity;

end

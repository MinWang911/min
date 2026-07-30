function [coded_bits, tx_sym, n_coded] = cc_encode_v4(tx_bits, cfg)
% CC_ENCODE_V4  Channel encode + optional 802.11a interleave + modulate
%
% INPUTS:
%   tx_bits — info bits (logical column vector)
%   cfg     — config struct:
%               cfg.scheme      : 'conv' (Rate 1/2, K=7) | 'ldpc' (DVB-S2, Rate 1/2)
%               cfg.modulation  : 'bpsk' | 'qpsk' | '16qam' | '64qam'
%               cfg.interleave  : (optional) true → apply 802.11a interleaver
%
% OUTPUTS:
%   coded_bits — encoded bits before interleaving (logical column vector)
%   tx_sym     — complex QAM symbols, unit average power
%   n_coded    — coded bit count before interleaving (pass to cc_decode_v4)

    scheme        = cfg.scheme;
    modulation    = cfg.modulation;
    do_interleave = isfield(cfg, 'interleave') && cfg.interleave;

    tx_bits = logical(tx_bits(:));
    M   = modulation_to_order(modulation);
    bps = log2(M);

    % Encoding
    switch lower(scheme)
        case 'conv'
            trellis    = poly2trellis(7, [171 133]);
            coded_bits = convenc(tx_bits, trellis);

        case 'ldpc'
            H      = dvbs2ldpc(1/2);
            cfgEnc = ldpcEncoderConfig(H);
            k_ldpc = cfgEnc.NumInformationBits;
            if length(tx_bits) < k_ldpc
                tx_bits = [tx_bits; false(k_ldpc - length(tx_bits), 1)];
            else
                tx_bits = tx_bits(1:k_ldpc);
            end
            coded_bits = ldpcEncode(tx_bits, cfgEnc);

        otherwise
            error('cc_encode_v4: cfg.scheme must be ''conv'' or ''ldpc''');
    end

    % n_coded recorded BEFORE interleaving — decoder needs this for LLR trimming
    n_coded = length(coded_bits);

    % 802.11a interleaving: spreads coded bits across subcarriers
    bits_to_map = coded_bits;
    if do_interleave
        bits_to_map = interleave_80211a(coded_bits, modulation);
    end

    % QAM mapping
    pad = mod(bps - mod(length(bits_to_map), bps), bps);
    coded_pad = [bits_to_map; false(pad, 1)];
    tx_sym = qammod(double(coded_pad), M, 'InputType', 'bit', 'UnitAveragePower', true);

end

function M = modulation_to_order(modulation)
    switch lower(modulation)
        case 'bpsk',  M = 2;
        case 'qpsk',  M = 4;
        case '16qam', M = 16;
        case '64qam', M = 64;
        otherwise, error('cc_encode_v4: unknown modulation ''%s''', modulation);
    end
end
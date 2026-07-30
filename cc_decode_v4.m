function rx_bits = cc_decode_v4(rx_eq, N0, cfg, n_coded)
% CC_DECODE_V4  Soft-demodulate + optional 802.11a deinterleave + decode
%
% INPUTS:
%   rx_eq   — equalised complex symbols (column vector)
%   N0      — noise variance per symbol (use mean(N0./abs(H_f).^2) after ZF)
%   cfg     — config struct:
%               cfg.scheme      : 'conv' | 'ldpc'
%               cfg.modulation  : 'bpsk' | 'qpsk' | '16qam' | '64qam'
%               cfg.interleave  : (optional) true → deinterleave LLRs before decoding
%   n_coded — (optional) from cc_encode_v4, for exact LLR trimming
%
% OUTPUT:
%   rx_bits — decoded info bits (logical column vector)

    if nargin < 4, n_coded = 0; end

    scheme        = cfg.scheme;
    modulation    = cfg.modulation;
    do_interleave = isfield(cfg, 'interleave') && cfg.interleave;

    rx_eq = rx_eq(:);
    M     = modulation_to_order(modulation);

    % Soft demodulation -> LLRs
    llr = qamdemod(rx_eq, M, 'OutputType', 'approxllr', ...
                   'UnitAveragePower', true, 'NoiseVariance', N0);
    llr = llr(:);

    % Deinterleave BEFORE trimming — full N_CBPS blocks must be intact
    if do_interleave
        llr = deinterleave_80211a(llr, modulation);
    end

    if n_coded > 0 && n_coded < length(llr)
        llr = llr(1:n_coded);
    end

    llr = max(min(llr, 20), -20);

    switch lower(scheme)
        case 'conv'
            trellis = poly2trellis(7, [171 133]);
            rx_bits = logical(vitdec(llr, trellis, 35, 'trunc', 'unquant'));

        case 'ldpc'
            H       = dvbs2ldpc(1/2);
            cfgDec  = ldpcDecoderConfig(H);
            rx_bits = logical(ldpcDecode(llr, cfgDec, 50));

        otherwise
            error('cc_decode_v4: cfg.scheme must be ''conv'' or ''ldpc''');
    end

end

function M = modulation_to_order(modulation)
    switch lower(modulation)
        case 'bpsk',  M = 2;
        case 'qpsk',  M = 4;
        case '16qam', M = 16;
        case '64qam', M = 64;
        otherwise, error('cc_decode_v4: unknown modulation ''%s''', modulation);
    end
end
function bits_out = interleave_80211n(bits_in, modulation)
% INTERLEAVE_80211N  IEEE 802.11n HT 20 MHz bit interleaver (single stream)
%
% 11n counterpart of interleave_80211a. The 11a interleaver used a 16-column
% layout tied to 48 data subcarriers (N_CBPS = 48*N_BPSC). 11n 20 MHz uses a
% 13-column layout tied to 52 data subcarriers:
%       N_COL = 13,  N_ROW = 4*N_BPSC,  N_CBPS = N_COL*N_ROW = 52*N_BPSC
% (IEEE 802.11n-2009, Sec 20.3.11.8.2). This spreads the coded bits across the
% 11n grid correctly - important for the frequency-diversity result.
%
% NOTE: the third (spatial-stream) permutation used when N_SS > 1 is NOT applied
% here. This is the PER-STREAM interleaver, which is what the 2x2 scripts need,
% since each spatial stream's bit stream is interleaved independently.
%
% INPUTS:
%   bits_in    - coded bits (logical column vector)
%   modulation - 'bpsk' | 'qpsk' | '16qam' | '64qam'
% OUTPUT:
%   bits_out   - interleaved bits (logical, length rounded UP to a multiple of
%                N_CBPS). Use n_coded from cc_encode_v4 to trim after
%                deinterleaving, exactly as with the 11a version.

    [N_CBPS, N_BPSC] = cbps_params_11n(modulation);
    perm = compute_perm_11n(N_CBPS, N_BPSC);

    bits_in = logical(bits_in(:));
    n       = length(bits_in);

    % Pad to a multiple of N_CBPS
    n_pad  = mod(N_CBPS - mod(n, N_CBPS), N_CBPS);
    padded = [bits_in; false(n_pad, 1)];
    n_sym  = length(padded) / N_CBPS;

    bits_out = false(size(padded));
    for sym = 1:n_sym
        idx     = (sym-1)*N_CBPS + (1:N_CBPS);
        out_blk = false(N_CBPS, 1);
        out_blk(perm) = padded(idx);   % scatter: input k -> output perm(k)
        bits_out(idx) = out_blk;
    end
    % Do NOT trim - returning full N_CBPS-aligned blocks keeps a partial last
    % block intact (same convention as interleave_80211a).
end

function [N_CBPS, N_BPSC] = cbps_params_11n(modulation)
    switch lower(modulation)
        case 'bpsk',  N_BPSC = 1;
        case 'qpsk',  N_BPSC = 2;
        case '16qam', N_BPSC = 4;
        case '64qam', N_BPSC = 6;
        otherwise, error('interleave_80211n: unknown modulation ''%s''', modulation);
    end
    N_CBPS = 52 * N_BPSC;          % 20 MHz: 52 data subcarriers
end

function perm = compute_perm_11n(N_CBPS, N_BPSC)
    N_COL = 13;                    % 20 MHz column count
    N_ROW = 4 * N_BPSC;            % = N_CBPS / N_COL
    s     = max(N_BPSC/2, 1);
    k = (0:N_CBPS-1)';
    % First permutation (column-row transpose)
    i = N_ROW * mod(k, N_COL) + floor(k / N_COL);
    % Second permutation (avoids long runs of same-reliability bits in high QAM)
    j = s * floor(i/s) + mod(i + N_CBPS - floor(N_COL * i / N_CBPS), s);
    perm = j + 1;                  % 1-based output index for each input bit k
end
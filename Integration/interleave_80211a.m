function bits_out = interleave_80211a(bits_in, modulation)
% INTERLEAVE_80211A  802.11a bit interleaver (IEEE 802.11a §17.3.5.6)
%
% INPUTS:
%   bits_in    — coded bits (logical column vector)
%   modulation — 'bpsk' | 'qpsk' | '16qam' | '64qam'
% OUTPUT:
%   bits_out   — interleaved bits (logical, length rounded UP to multiple of N_CBPS)
%
% NOTE: output length >= input length (padded to N_CBPS boundary).
%       Caller should use n_coded from cc_encode_v4 to trim after deinterleaving.

    [N_CBPS, N_BPSC] = cbps_params(modulation);
    s    = max(N_BPSC / 2, 1);
    perm = compute_perm(N_CBPS, s);

    bits_in = logical(bits_in(:));
    n       = length(bits_in);

    % Pad to multiple of N_CBPS
    n_pad  = mod(N_CBPS - mod(n, N_CBPS), N_CBPS);
    padded = [bits_in; false(n_pad, 1)];
    n_sym  = length(padded) / N_CBPS;

    bits_out = false(size(padded));
    for sym = 1:n_sym
        idx     = (sym-1)*N_CBPS + (1:N_CBPS);
        out_blk = false(N_CBPS, 1);
        out_blk(perm) = padded(idx);   % scatter: input k → output perm(k)
        bits_out(idx) = out_blk;
    end

    % Do NOT trim — returning full N_CBPS-aligned block keeps partial
    % last block intact. Trimming here would corrupt partial-block bits
    % because the permutation mixes real bits with padding zeros.

end

function [N_CBPS, N_BPSC] = cbps_params(modulation)
    switch lower(modulation)
        case 'bpsk',  N_CBPS = 48;  N_BPSC = 1;
        case 'qpsk',  N_CBPS = 96;  N_BPSC = 2;
        case '16qam', N_CBPS = 192; N_BPSC = 4;
        case '64qam', N_CBPS = 288; N_BPSC = 6;
        otherwise, error('interleave_80211a: unknown modulation ''%s''', modulation);
    end
end

function perm = compute_perm(N_CBPS, s)
    k    = (0:N_CBPS-1)';
    i    = (N_CBPS/16) * mod(k, 16) + floor(k/16);
    j    = s * floor(i/s) + mod(i + N_CBPS - floor(16*i/N_CBPS), s);
    perm = j + 1;
end
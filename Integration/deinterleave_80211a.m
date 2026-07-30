function llr_out = deinterleave_80211a(llr_in, modulation)
% DEINTERLEAVE_80211A  Inverse of interleave_80211a — applied to LLRs at RX
%
% INPUTS:
%   llr_in     — interleaved LLRs (real column vector, from qamdemod)
%   modulation — 'bpsk' | 'qpsk' | '16qam' | '64qam'
% OUTPUT:
%   llr_out    — deinterleaved LLRs (real, same length as llr_in)

    [N_CBPS, N_BPSC] = cbps_params(modulation);
    s    = max(N_BPSC / 2, 1);
    perm = compute_perm(N_CBPS, s);   % same permutation as interleaver

    llr_in = double(llr_in(:));
    n      = length(llr_in);

    n_pad  = mod(N_CBPS - mod(n, N_CBPS), N_CBPS);
    padded = [llr_in; zeros(n_pad, 1)];
    n_sym  = length(padded) / N_CBPS;

    llr_out = zeros(size(padded));
    for sym = 1:n_sym
        idx          = (sym-1)*N_CBPS + (1:N_CBPS);
        llr_out(idx) = padded(idx(perm));   % gather: read perm(k) → output k
    end

    llr_out = llr_out(1:n);
end

function [N_CBPS, N_BPSC] = cbps_params(modulation)
    switch lower(modulation)
        case 'bpsk',  N_CBPS = 48;  N_BPSC = 1;
        case 'qpsk',  N_CBPS = 96;  N_BPSC = 2;
        case '16qam', N_CBPS = 192; N_BPSC = 4;
        case '64qam', N_CBPS = 288; N_BPSC = 6;
        otherwise, error('deinterleave_80211a: unknown modulation ''%s''', modulation);
    end
end

function perm = compute_perm(N_CBPS, s)
    k    = (0:N_CBPS-1)';
    i    = (N_CBPS/16) * mod(k, 16) + floor(k/16);
    j    = s * floor(i/s) + mod(i + N_CBPS - floor(16*i/N_CBPS), s);
    perm = j + 1;
end
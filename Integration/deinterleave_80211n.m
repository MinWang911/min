function llr_out = deinterleave_80211n(llr_in, modulation)
% DEINTERLEAVE_80211N  Inverse of interleave_80211n - applied to LLRs at RX.
%
% IEEE 802.11n HT 20 MHz, single spatial stream.
%   N_COL = 13,  N_ROW = 4*N_BPSC,  N_CBPS = 52*N_BPSC
% Uses the SAME permutation as interleave_80211n and gathers it back, exactly
% as deinterleave_80211a does for the 11a grid.
%
% INPUTS:
%   llr_in     - interleaved LLRs (real column vector, from qamdemod)
%   modulation - 'bpsk' | 'qpsk' | '16qam' | '64qam'
% OUTPUT:
%   llr_out    - deinterleaved LLRs (real, same length as llr_in)

    [N_CBPS, N_BPSC] = cbps_params_11n(modulation);
    perm = compute_perm_11n(N_CBPS, N_BPSC);

    llr_in = double(llr_in(:));
    n      = length(llr_in);

    n_pad  = mod(N_CBPS - mod(n, N_CBPS), N_CBPS);
    padded = [llr_in; zeros(n_pad, 1)];
    n_sym  = length(padded) / N_CBPS;

    llr_out = zeros(size(padded));
    for sym = 1:n_sym
        idx          = (sym-1)*N_CBPS + (1:N_CBPS);
        llr_out(idx) = padded(idx(perm));   % gather: read perm(k) -> output k
    end
    llr_out = llr_out(1:n);
end

function [N_CBPS, N_BPSC] = cbps_params_11n(modulation)
    switch lower(modulation)
        case 'bpsk',  N_BPSC = 1;
        case 'qpsk',  N_BPSC = 2;
        case '16qam', N_BPSC = 4;
        case '64qam', N_BPSC = 6;
        otherwise, error('deinterleave_80211n: unknown modulation ''%s''', modulation);
    end
    N_CBPS = 52 * N_BPSC;          % 20 MHz: 52 data subcarriers
end

function perm = compute_perm_11n(N_CBPS, N_BPSC)
    N_COL = 13;                    % 20 MHz column count
    N_ROW = 4 * N_BPSC;            % = N_CBPS / N_COL
    s     = max(N_BPSC/2, 1);
    k = (0:N_CBPS-1)';
    i = N_ROW * mod(k, N_COL) + floor(k / N_COL);
    j = s * floor(i/s) + mod(i + N_CBPS - floor(N_COL * i / N_CBPS), s);
    perm = j + 1;
end
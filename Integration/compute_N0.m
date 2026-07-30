function N0 = compute_N0(SNR_dB, code_rate, M)
% COMPUTE_N0  Noise power per symbol for M-QAM channel coding
%
% INPUTS:
%   SNR_dB    — Eb/N0 in dB
%   code_rate — e.g. 0.5 for Rate 1/2
%   M         — modulation order: 2, 4, 16, or 64
%
% OUTPUT:
%   N0  — noise power per symbol (pass to cc_decode_v3)
%         After ZF equalisation: pass mean(N0 ./ abs(H_f).^2) instead

    N0 = 1 / (10^(SNR_dB/10) * log2(M) * code_rate);

end

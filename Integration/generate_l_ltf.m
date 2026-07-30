function [L_LTF_freq, ltf_seq] = generate_l_ltf(cfg)
% GENERATE_L_LTF Generates the standard IEEE 802.11a/p L-LTF sequence.
%
% Input:
%   cfg          - Structure containing OFDM parameters (N_fft, etc.)
%
% Outputs:
%   L_LTF_freq   - 64x1 (or N_fft) frequency-domain vector with L-LTF mapped
%   ltf_seq      - The raw 52-length +1/-1 standard L-LTF sequence (Table 17-6)

    % Standard 802.11a/p L-LTF sequence for subcarriers -26 to 26 (including DC at index 0)
    % Note: 0 at DC (index 27 in shifted grid)
    ltf_seq = [1, 1, -1, -1, 1, 1, -1, 1, -1, 1, 1, 1, 1, 1, 1, -1, -1, 1, 1, -1, 1, -1, 1, 1, 1, 1, ...
               0, ... % DC Subcarrier
               1, -1, -1, 1, 1, -1, 1, -1, 1, -1, -1, -1, -1, -1, 1, 1, -1, -1, 1, -1, 1, -1, 1, 1, 1, 1];

    % Initialize full FFT grid
    L_LTF_freq = zeros(cfg.N_fft, 1);
    
    % Map L-LTF sequence to standard subcarriers [-26:26]
    % In MATLAB index: -26 to 26 maps to indices 7 to 59 (using N_fft=64)
    subcarrier_indices_80211 = -26:26;
    to_matlab_idx = @(x) x + (cfg.N_fft/2) + 1;
    mapped_indices = to_matlab_idx(subcarrier_indices_80211);
    
    L_LTF_freq(mapped_indices) = ltf_seq.';
end
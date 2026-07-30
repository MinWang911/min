function LTF_shifted = generate_known_ltf(cfg)
% GENERATE_KNOWN_LTF
% Generate known 802.11a LTF sequence in fftshifted 64-subcarrier form.
%
% Input:
%   cfg : system configuration structure
%
% Required cfg fields:
%   cfg.N_fft
%
% Output:
%   LTF_shifted : known LTF sequence on full fftshifted OFDM grid

    if ~isfield(cfg, 'N_fft')
        error('Missing cfg field: N_fft');
    end

    Nfft = cfg.N_fft;

    if Nfft ~= 64
        error('This LTF sequence is for 64-point 802.11a OFDM.');
    end

    % LTF values for subcarriers -26 to +26.
    % DC subcarrier is 0.
    ltf_53 = [ ...
         1;  1; -1; -1;  1;  1; -1;  1; -1;  1;  1;  1;  1; ...
         1;  1; -1; -1;  1;  1; -1;  1; -1;  1;  1;  1;  1; ...
         0; ...
         1; -1; -1;  1;  1; -1;  1; -1;  1; -1; -1; -1; -1; ...
        -1;  1;  1; -1; -1;  1; -1;  1; -1;  1;  1;  1;  1];

    sc = (-26:26).';

    LTF_shifted = zeros(Nfft, 1);

    scToIndex = @(x) x + Nfft/2 + 1;

    LTF_shifted(scToIndex(sc)) = ltf_53;

end

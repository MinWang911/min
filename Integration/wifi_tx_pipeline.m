function [tx_signal, tx_bits] = wifi_tx_pipeline(cfg, tx_bits)
% WIFI_TX_PIPELINE IEEE 802.11a/p Transmitter with L-LTF Preamble Insertion
%
% Inputs:
%   cfg       - Structure containing OFDM parameters (N_fft, N_cp, ModOrder, etc.)
%   tx_bits   - (Optional) Channel encoded bitstream (logical/double/uint8).
%               If omitted, auto-generates random bits.
%
% Outputs:
%   tx_signal - Time-domain vector with Preamble (2 L-LTFs) + Data symbols
%   tx_bits   - The actual bits transmitted (returned as 'double')

    %% --- Step 1: Bitstream Handling & Dynamic Allocation ---
    bitsPerSymbol = cfg.N_data * log2(cfg.ModOrder); 
    expectedLength = bitsPerSymbol * cfg.numSymbols;

    if nargin < 2 || isempty(tx_bits)
        tx_bits = randi([0 1], expectedLength, 1);
    else
        tx_bits = tx_bits(:);
        if length(tx_bits) ~= expectedLength
            error('WIFI_TX:BitLengthMismatch', ...
                'Input bitstream length (%d) does not match expected capacity (%d bits).', ...
                length(tx_bits), expectedLength);
        end
    end
    tx_bits = double(tx_bits); 

    %% --- Step 2: Generate and Process L-LTF Preamble ---
    % Call the newly modularized L-LTF sequence generator
    [L_LTF_freq, ~] = generate_l_ltf(cfg);
    
    % Power normalization for L-LTF (Standard 52 active subcarriers)
    norm_factor_preamble = sqrt(cfg.N_fft / 52); 
    L_LTF_freq_norm = L_LTF_freq * norm_factor_preamble;
    
    % Convert 1 L-LTF symbol to Time-Domain via IFFT
    ltf_time_no_cp = ifft(ifftshift(L_LTF_freq_norm, 1), cfg.N_fft);
    
    % Add Cyclic Prefix for L-LTF
    ltf_cp = ltf_time_no_cp(end-cfg.N_cp+1:end);
    ltf_symbol_time = [ltf_cp; ltf_time_no_cp]; % Length: N_fft + N_cp (80 samples)
    
    % Duplicate to create 2 consecutive L-LTF symbols as required by the standard
    preamble_time = [ltf_symbol_time; ltf_symbol_time]; % Total: 160 samples

    %% --- Step 3: Standard Constellation Modulation for Payload ---
    tx_modulated = qammod(tx_bits, cfg.ModOrder, 'InputType', 'bit', 'UnitAveragePower', true);

    %% --- Step 4: IEEE 802.11a/p Pilot PRBS Generator ---
    state = ones(1, 7); 
    p_seq = zeros(1, cfg.numSymbols);
    for n = 1:cfg.numSymbols
        p_seq(n) = state(7); 
        feedback = xor(state(7), state(4)); 
        state = [feedback, state(1:6)];
    end
    pilot_polarity = 1 - 2 * p_seq; 
    pilot_base = [1; 1; 1; -1]; 

    %% --- Step 5: Payload Subcarrier Mapping & Power Normalization ---
    X_freq = zeros(cfg.N_fft, cfg.numSymbols);
    X_freq(cfg.data_idx, :)  = reshape(tx_modulated, cfg.N_data, cfg.numSymbols);
    X_freq(cfg.pilot_idx, :) = pilot_base * pilot_polarity;

    norm_factor_payload = sqrt(cfg.N_fft / (cfg.N_data + cfg.N_pilot)); 
    X_freq = X_freq * norm_factor_payload;

    %% --- Step 6: Payload IFFT & CP Insertion ---
    payload_time_no_cp = ifft(ifftshift(X_freq, 1), cfg.N_fft);
    payload_cp = payload_time_no_cp(end-cfg.N_cp+1:end, :);
    payload_matrix = [payload_cp; payload_time_no_cp]; 
    payload_time = payload_matrix(:); % Serialize payload symbols

    %% --- Step 7: Frame Assembly (Prepend Preamble to Payload) ---
    % Assemble: [ 2x L-LTF Symbols (160 samples) ] + [ Payload Symbols ]
    tx_signal = [preamble_time; payload_time];
end

function Hhat = mimo_estimate_channel(rxTrain, C, P)
%MIMO_ESTIMATE_CHANNEL  Least-squares MIMO channel estimate (grid domain).
%
% Real (non-genie) frequency-domain MIMO channel estimation. Each Tx antenna
% transmits the known pilot P(sc) weighted by an orthogonal cover code C across
% Ntrain training OFDM symbols. Because the channel is held constant over the
% frame, the receiver separates the Ntx antenna contributions per subcarrier by
% inverting C in a least-squares sense.
%
% This is the MIMO analogue of the SISO LTF estimator (estimate_channel.m),
% which assumes a single Tx antenna and therefore cannot be reused directly.
%
% INPUTS
%   rxTrain : [Nsc x Ntrain x Nrx]  received training grid (after the channel)
%   C       : [Ntrain x Ntx]        orthogonal cover, e.g. hadamard(Ntx); Ntrain>=Ntx
%   P       : [Nsc x 1]             known pilot value per subcarrier (0 on nulls)
%
% OUTPUT
%   Hhat    : [Nrx x Ntx x Nsc]     LS channel estimate (same layout as ofdm_channel H)
%
% Model, per subcarrier sc and Rx antenna r:
%   y_r = P(sc) * C * h_r + noise,   h_r = squeeze(H(r,:,sc)).'
%   => h_r = pinv(C) * (y_r / P(sc))
%
% With a Hadamard cover and Ntrain = Ntx this is the standard orthogonal-pilot
% MIMO LS estimate (C'C = Ntx*I).

    [Nsc, Ntrain, Nrx] = size(rxTrain);
    Ntx = size(C, 2);
    assert(size(C,1) == Ntrain, 'C must have Ntrain rows.');
    assert(numel(P)  == Nsc,    'P must have Nsc entries.');

    Cpinv = pinv(C);                          % [Ntx x Ntrain]
    Hhat  = complex(zeros(Nrx, Ntx, Nsc));

    for sc = 1:Nsc
        if P(sc) == 0, continue; end
        Y   = reshape(rxTrain(sc,:,:), Ntrain, Nrx);   % [Ntrain x Nrx]
        Hsc = Cpinv * (Y / P(sc));                     % [Ntx x Nrx], col r = h_r
        Hhat(:,:,sc) = Hsc.';                          % [Nrx x Ntx]
    end
end

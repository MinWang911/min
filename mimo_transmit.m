function [txAntennaGrid, meta] = mimo_transmit(txLayerGrid, cfg)
%MIMO_TRANSMIT Map data layers to transmit antennas.
%
% Spatial multiplexing input: [Nsc x Nsym x Nstream]
% Alamouti input:              [Nsc x Nblock x 2] (s1,s2 pairs)
% Output:                      [Nsc x Nslot x Ntx]

    required = {'mode','numTx','numRx','numStreams','precoder'};
    assert(all(isfield(cfg, required)), 'MIMO configuration is incomplete.');

    switch lower(cfg.mode)
        case 'spatial_multiplexing'
            [numSC, numSym, numStreams] = size(txLayerGrid);
            assert(numStreams == cfg.numStreams, ...
                'Layer count must equal cfg.numStreams.');
            validate_precoder(cfg.precoder, cfg, numSC, numSym);

            txAntennaGrid = complex(zeros(numSC, numSym, cfg.numTx));
            for sc = 1:numSC
                for sym = 1:numSym
                    W = precoder_at(cfg.precoder, sc, sym);
                    x = W * reshape(txLayerGrid(sc,sym,:), [], 1);
                    txAntennaGrid(sc,sym,:) = reshape(x, 1, 1, []);
                end
            end
            meta.numInputSlots = numSym;
            meta.numOutputSlots = numSym;
            meta.codeRate = 1;

        case 'alamouti'
            [numSC, numBlocks, numLayers] = size(txLayerGrid);
            assert(cfg.numTx == 2 && numLayers == 2, ...
                'Alamouti input must contain two-symbol pairs for two transmit antennas.');
            s1 = txLayerGrid(:,:,1);
            s2 = txLayerGrid(:,:,2);
            txAntennaGrid = complex(zeros(numSC, 2*numBlocks, 2));
            txAntennaGrid(:,1:2:end,1) =  s1 / sqrt(2);
            txAntennaGrid(:,1:2:end,2) =  s2 / sqrt(2);
            txAntennaGrid(:,2:2:end,1) = -conj(s2) / sqrt(2);
            txAntennaGrid(:,2:2:end,2) =  conj(s1) / sqrt(2);
            meta.numInputSlots = numBlocks;
            meta.numOutputSlots = 2*numBlocks;
            meta.codeRate = 1;

        otherwise
            error('Unsupported MIMO mode: %s', cfg.mode);
    end
    meta.mode = cfg.mode;
    meta.numTx = cfg.numTx;
end

function validate_precoder(P, cfg, numSC, numSym)
    assert(size(P,1) == cfg.numTx && size(P,2) == cfg.numStreams, ...
        'Precoder first dimensions must be [numTx x numStreams].');
    assert(size(P,3) == 1 || size(P,3) == numSC, ...
        'Precoder subcarrier dimension must be 1 or Nsc.');
    assert(size(P,4) == 1 || size(P,4) == numSym, ...
        'Precoder symbol dimension must be 1 or Nsym.');
    assert(all(isfinite(P(:))), 'Precoder must contain finite values.');
end

function W = precoder_at(P, sc, sym)
    W = P(:,:,min(sc,size(P,3)),min(sym,size(P,4)));
end

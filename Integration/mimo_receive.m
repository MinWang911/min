function [rxLayerGrid, meta] = mimo_receive(rxGrid, Hhat, noiseVariance, cfg)
%MIMO_RECEIVE Recover spatial layers using perfect or estimated CSI.
%
% rxGrid: [Nsc x Nslot x Nrx]
% Hhat:   [Nrx x Ntx x Nsc] or [Nrx x Ntx x Nsc x Nslot]

    validate_inputs(rxGrid, Hhat, cfg);
    [numSC, numSlots, ~] = size(rxGrid);

    switch lower(cfg.mode)
        case 'spatial_multiplexing'
            rxLayerGrid = complex(zeros(numSC, numSlots, cfg.numStreams));
            for sc = 1:numSC
                for sym = 1:numSlots
                    H = channel_at(Hhat, sc, sym) * precoder_at(cfg.precoder, sc, sym);
                    y = reshape(rxGrid(sc,sym,:), [], 1);
                    nv = noise_at(noiseVariance, sc, sym);
                    xHat = mimo_detect(y, H, nv, cfg.detector, cfg.constellation);
                    rxLayerGrid(sc,sym,:) = reshape(xHat,1,1,[]);
                end
            end
            meta.numDetectedLayers = cfg.numStreams;

        case 'alamouti'
            assert(mod(numSlots,2) == 0, 'Alamouti requires an even number of slots.');
            rxLayerGrid = complex(zeros(numSC, numSlots/2, 2));
            tol = cfg.alamoutiTolerance;
            for sc = 1:numSC
                for slot = 1:2:numSlots
                    H1 = channel_at(Hhat, sc, slot);
                    H2 = channel_at(Hhat, sc, slot+1);
                    assert(norm(H1-H2,'fro') <= tol*max(1,norm(H1,'fro')), ...
                        'Channel must remain constant across each Alamouti slot pair.');
                    h1 = H1(:,1)/sqrt(2);
                    h2 = H1(:,2)/sqrt(2);
                    y1 = reshape(rxGrid(sc,slot,:), [], 1);
                    y2 = reshape(rxGrid(sc,slot+1,:), [], 1);
                    den = sum(abs(h1).^2 + abs(h2).^2) + eps;
                    s1 = sum(conj(h1).*y1 + h2.*conj(y2))/den;
                    s2 = sum(conj(h2).*y1 - h1.*conj(y2))/den;
                    rxLayerGrid(sc,(slot+1)/2,:) = reshape([s1;s2],1,1,2);
                end
            end
            meta.numDetectedLayers = 2;

        otherwise
            error('Unsupported MIMO mode: %s', cfg.mode);
    end
    meta.mode = cfg.mode;
    meta.detector = cfg.detector;
end

function validate_inputs(rxGrid, H, cfg)
    required = {'mode','numTx','numRx','numStreams','precoder','detector','constellation'};
    assert(all(isfield(cfg,required)), 'MIMO configuration is incomplete.');
    assert(size(rxGrid,3) == cfg.numRx, 'rxGrid receiver dimension is incorrect.');
    assert(size(H,1) == cfg.numRx && size(H,2) == cfg.numTx, ...
        'Hhat antenna dimensions do not match cfg.');
    assert(size(H,3) == size(rxGrid,1), 'Hhat subcarrier dimension is incorrect.');
    assert(size(H,4) == 1 || size(H,4) == size(rxGrid,2), ...
        'Hhat slot dimension must be 1 or match rxGrid.');
end

function Hsc = channel_at(H, sc, sym)
    Hsc = H(:,:,sc,min(sym,size(H,4)));
end

function W = precoder_at(P, sc, sym)
    W = P(:,:,min(sc,size(P,3)),min(sym,size(P,4)));
end

function nv = noise_at(noiseVariance, sc, sym)
    if isscalar(noiseVariance)
        nv = noiseVariance;
    else
        nv = noiseVariance(sc,sym);
    end
end

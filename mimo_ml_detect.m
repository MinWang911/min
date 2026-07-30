function xHat = mimo_ml_detect(y, H, constellation)
%MIMO_ML_DETECT Exhaustive ML detector for small MIMO constellations.
% Complexity is M^Ns, so this is intended for small Ns and low-order QAM.

    constellation = constellation(:).';
    numStreams = size(H,2);
    M = numel(constellation);
    candidateCount = M^numStreams;
    candidates = complex(zeros(numStreams, candidateCount));
    indices = 0:candidateCount-1;
    for stream = 1:numStreams
        digit = mod(floor(indices / M^(stream-1)), M) + 1;
        candidates(stream,:) = constellation(digit);
    end
    residual = y(:) - H*candidates;
    [~, best] = min(sum(abs(residual).^2, 1));
    xHat = candidates(:,best);
end

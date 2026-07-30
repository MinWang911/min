function metrics = mimo_evaluate(txBits,rxBits,txSymbols,rxSymbols,numResourceElements)
%MIMO_EVALUATE Calculate BER, RMS EVM and net bits per resource element.

    txBits = txBits(:); rxBits = rxBits(:);
    txSymbols = txSymbols(:); rxSymbols = rxSymbols(:);
    bitCount = min(numel(txBits),numel(rxBits));
    symbolCount = min(numel(txSymbols),numel(rxSymbols));
    assert(bitCount > 0 && symbolCount > 0, 'Evaluation inputs must not be empty.');
    metrics.numBits = bitCount;
    metrics.bitErrors = sum(txBits(1:bitCount) ~= rxBits(1:bitCount));
    metrics.ber = metrics.bitErrors/bitCount;
    metrics.rmsEVM = sqrt(sum(abs(rxSymbols(1:symbolCount)-txSymbols(1:symbolCount)).^2) / ...
        max(sum(abs(txSymbols(1:symbolCount)).^2),eps));
    if nargin >= 5 && ~isempty(numResourceElements)
        metrics.netBitsPerRE = bitCount/numResourceElements;
    else
        metrics.netBitsPerRE = NaN;
    end
end

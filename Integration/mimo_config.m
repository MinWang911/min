function cfg = mimo_config(mode, numTx, numRx, numStreams)
%MIMO_CONFIG Create a unified configuration for the integrated MIMO module.
%
% cfg = mimo_config('spatial_multiplexing', 2, 2, 2)
% cfg = mimo_config('alamouti', 2, 2)

    if nargin < 1 || isempty(mode), mode = 'spatial_multiplexing'; end
    if nargin < 2 || isempty(numTx), numTx = 2; end
    if nargin < 3 || isempty(numRx), numRx = 2; end

    mode = lower(strrep(mode, ' ', '_'));
    if strcmp(mode, 'smux'), mode = 'spatial_multiplexing'; end
    if strcmp(mode, 'diversity'), mode = 'alamouti'; end

    validateattributes(numTx, {'numeric'}, {'scalar','integer','positive'});
    validateattributes(numRx, {'numeric'}, {'scalar','integer','positive'});

    switch mode
        case 'spatial_multiplexing'
            if nargin < 4 || isempty(numStreams)
                numStreams = min(numTx, numRx);
            end
            assert(numStreams <= min(numTx, numRx), ...
                'Spatial multiplexing needs numStreams <= min(numTx,numRx).');
        case 'alamouti'
            assert(numTx == 2, 'Alamouti mode requires exactly two transmit antennas.');
            numStreams = 2;  % s1 and s2 are carried over two time slots
        otherwise
            error('Unsupported MIMO mode: %s', mode);
    end

    cfg.mode = mode;
    cfg.numTx = numTx;
    cfg.numRx = numRx;
    cfg.numStreams = numStreams;
    cfg.detector = 'mmse';
    cfg.precoder = eye(numTx, numStreams);
    cfg.constellation = [1+1i; 1-1i; -1+1i; -1-1i] / sqrt(2);
    cfg.alamoutiTolerance = 1e-8;
end

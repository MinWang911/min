function xHat = mimo_detect(y, H, noiseVariance, detector, constellation)
%MIMO_DETECT Detect one spatial-multiplexing resource element.
% detector: 'zf', 'mmse', or 'ml'.

    if nargin < 4 || isempty(detector), detector = 'mmse'; end
    if nargin < 5, constellation = []; end
    y = y(:);
    numStreams = size(H,2);

    switch lower(detector)
        case 'zf'
            xHat = pinv(H) * y;
        case 'mmse'
            validateattributes(noiseVariance, {'numeric'}, {'scalar','real','nonnegative'});
            xHat = (H'*H + noiseVariance*eye(numStreams)) \ (H'*y);
        case 'ml'
            assert(~isempty(constellation), 'ML detection requires cfg.constellation.');
            xHat = mimo_ml_detect(y, H, constellation);
        otherwise
            error('Unsupported detector "%s". Use ZF, MMSE, or ML.', detector);
    end
end

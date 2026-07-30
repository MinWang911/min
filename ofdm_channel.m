function [rxSignal, H, info] = ofdm_channel(txSignal, cfg)
%OFDM_CHANNEL  Channel model for the integrated system.
%              Channel Modelling (Role 4). Wahab. Group G3005.
%
% ============================================================================
% THE SINGLE ENTRY POINT for the channel. Call this from the integration
% script. You do not need to read my analysis scripts.
%
%   [rxSignal, H, info] = ofdm_channel(txSignal, cfg)
%
% Handles SISO and MIMO, and two operating domains. Returns the faded + noisy
% signal, the TRUE channel response H, and a diagnostics struct.
%
% H is ground truth. Use it as perfect CSI for a baseline, or compare your own
% estimate against it to measure estimation error. Returning it is not a
% shortcut, it is how the estimator gets scored.
%
% ---------------------------------------------------------------------------
% TWO DOMAINS  (cfg.domain)
%
%   'time'   USE THIS FOR INTEGRATION
%       txSignal : [nSamples x Ntx] time-domain, CP already attached.
%       - LINEAR convolution, so real ISI appears when the maximum delay
%         exceeds the CP (info.cpExceeded).
%       - Per-tap time-varying gains (sum-of-sinusoids Jakes), so Doppler
%         produces genuine intra-symbol variation and therefore real ICI.
%       - The Rician LOS ray carries its own deterministic Doppler shift, so
%         a strong LOS makes ICI WORSE for a fast terminal (e.g. LEO), not
%         better. Freezing the LOS understates the satellite floor by ~2x.
%       - Noise added at the RECEIVER FRONT END, before any FFT, so it
%         corrupts PILOTS and DATA consistently. Required to test a
%         pilot-based estimator fairly.
%
%   'grid'   fast, for isolated characterisation with perfect equalisation
%       txSignal : [nSc x nSym x Ntx] frequency-domain grid.
%       - Applies H per subcarrier (circular convolution assumption).
%       - CANNOT produce ISI (delays wrap) and CANNOT produce ICI (the
%         channel is constant within a symbol by construction).
%       - Equivalent to 'time' only when max delay < CP and CSI is perfect.
%
% ---------------------------------------------------------------------------
% NOISE CONVENTION - the single easiest thing to get wrong
%
%   noiseVar = (measured TRANSMITTED power * Ntx) / snrLin
%
%   Measure the TRANSMITTED signal. Two wrong alternatives, both tempting:
%
%   (a) Assuming unit power. After the IFFT with 48/64 subcarriers loaded the
%       time-domain power is ~0.0117, not 1.0. Assuming 1.0 understates the
%       delivered SNR by ~19.3 dB.
%
%   (b) Measuring the FADED (received) power. This rescales the noise to each
%       fade, so a deep fade gets proportionally less noise and the fading
%       penalty cancels itself out. BER then tracks AWGN, not Rayleigh. This
%       is the classic awgn(...,'measured') trap. Verified: at Eb/N0 = 0 dB it
%       returns 7.8e-2 (AWGN theory) instead of 1.46e-1 (Rayleigh theory).
%
%   The x Ntx factor keeps the per-receive-antenna SNR convention: each Rx
%   antenna sums Ntx streams and E|H|^2 = 1, so mean received power = Ntx*txPow.
%   This matches mimo_channel's noiseVar = Ntx/snrLin when txPow = 1.
%
% ---------------------------------------------------------------------------
% SNR CONVENTION  (cfg.snrMode) - set explicitly, do not guess
%   'snr'   (default) cfg.snrdB is received SNR per receive antenna.
%           Use this for MIMO, matching mimo_channel.
%   'ebno'  cfg.snrdB is Eb/N0; converted using bitsPerSym and the data
%           subcarrier occupancy. Use this to compare against my Stage 1-5
%           single-antenna curves and against flat-Rayleigh theory.
%
% ---------------------------------------------------------------------------
% CONFIG (cfg)
%   REQUIRED
%     .domain            'time' | 'grid'
%     .snrdB             see snrMode above
%     .sampleRate        Hz (e.g. 20e6)
%
%   PICK AN ENVIRONMENT (easiest)
%     .environment       'indoor_home' | 'indoor_office' | 'hotspot' |
%                        'urban' | 'highway' | 'satellite'
%     ... or set .rmsDelaySpread_s, .K, .speed_kmh directly (these override).
%
%   OFDM
%     .nFFT     default 64        .nCP        default 16
%     .nDataSC  default 48        .bitsPerSym default 2   (both: 'ebno' only)
%     .fc       default 5.2e9
%     .scIdx    subcarrier indices carried in a 'grid' txSignal.
%               default 1:nSc. For 802.11a pass the real data indices so H is
%               evaluated at the right frequency bins.
%
%   MIMO (omit for SISO)
%     .mimo.numTx .mimo.numRx
%     .mimo.rho          spatial correlation, or .rhoTx / .rhoRx
%
%   OPTIONS
%     .holdOverSymbols   default true. Freezes the channel across the frame.
%                        REQUIRED for Alamouti (mimo_receive asserts H is
%                        constant over each 2-symbol block). Set false to let
%                        Doppler vary the channel symbol-to-symbol; H is then
%                        returned per symbol (4-D for MIMO).
%     .addNoise default true   .seed rng seed   .verbose default false
%
% ---------------------------------------------------------------------------
% OUTPUT
%   rxSignal : 'time' -> [nSamples x Nrx]     'grid' -> [nSc x nSym x Ntx->Nrx]
%   H        : SISO 'time' -> [nFFT x nSym]   (per-symbol average channel)
%              SISO 'grid' -> [nSc  x nSym]
%              MIMO 'time' -> [Nrx x Ntx x nFFT x nSym]
%              MIMO 'grid' -> [Nrx x Ntx x nSc  x nSym]
%              If holdOverSymbols (default) the symbol dimension is dropped,
%              giving the 3-D form mimo_receive expects.
%   info     : noiseVariance, measuredTxPow, cpExceeded, powerBeyondCP,
%              isiSeverity, bindingLimit,
%              coherenceBW_Hz, BcOverSC, fd_Hz, coherenceTime_s, TcOverTsym,
%              delayTaps_s, pdp, rhoTx, rhoRx, domain, snrMode
% ============================================================================

% ---------------- defaults ----------------
cfg = setdef(cfg,'domain','time');
cfg = setdef(cfg,'snrMode','snr');
cfg = setdef(cfg,'nFFT',64);
cfg = setdef(cfg,'nCP',16);
cfg = setdef(cfg,'nDataSC',48);
cfg = setdef(cfg,'bitsPerSym',2);
cfg = setdef(cfg,'fc',5.2e9);
cfg = setdef(cfg,'addNoise',true);
cfg = setdef(cfg,'verbose',false);
cfg = setdef(cfg,'holdOverSymbols',true);
assert(isfield(cfg,'snrdB'),      'ofdm_channel: cfg.snrdB is required.');
assert(isfield(cfg,'sampleRate'), 'ofdm_channel: cfg.sampleRate is required.');
if isfield(cfg,'seed') && ~isempty(cfg.seed), rng(cfg.seed); end

% ---------------- environment presets (Stage 5) ----------------
if isfield(cfg,'environment') && ~isempty(cfg.environment)
    p = envPreset(cfg.environment);
    cfg = setdef(cfg,'rmsDelaySpread_s',p.rms);
    cfg = setdef(cfg,'K',p.K);
    cfg = setdef(cfg,'speed_kmh',p.v);
end
cfg = setdef(cfg,'rmsDelaySpread_s',30e-9);
cfg = setdef(cfg,'K',0);
cfg = setdef(cfg,'speed_kmh',0);

% ---------------- MIMO or SISO ----------------
isMimo = isfield(cfg,'mimo') && ~isempty(cfg.mimo);
if isMimo
    Ntx = cfg.mimo.numTx;  Nrx = cfg.mimo.numRx;
    rhoTx = 0; rhoRx = 0;
    if isfield(cfg.mimo,'rho'),   rhoTx = cfg.mimo.rho; rhoRx = cfg.mimo.rho; end
    if isfield(cfg.mimo,'rhoTx'), rhoTx = cfg.mimo.rhoTx; end
    if isfield(cfg.mimo,'rhoRx'), rhoRx = cfg.mimo.rhoRx; end
else
    Ntx = 1; Nrx = 1; rhoTx = 0; rhoRx = 0;
end

% ---------------- derived physics ----------------
c     = 3e8;
sr    = cfg.sampleRate;
rms_s = cfg.rmsDelaySpread_s;
Tsym  = (cfg.nFFT + cfg.nCP)/sr;
cpDur = cfg.nCP/sr;
fd    = (cfg.speed_kmh/3.6)*cfg.fc/c;

% ---------------- delay profile (exponential PDP) ----------------
nTaps     = 6;
tau       = linspace(0, 5*rms_s, nTaps);
pdp       = exp(-tau/rms_s); pdp = pdp/sum(pdp);
delaySamp = round(tau*sr);
maxDelay  = max(tau);

% ISI depends on the MAXIMUM excess delay, not the RMS. But a boolean flag
% over-simplifies: what actually sets the ISI floor is HOW MUCH CHANNEL POWER
% falls outside the guard interval. Report both.
cpExceeded    = maxDelay > cpDur;
% A tap landing exactly ON the CP boundary is INSIDE the guard interval.
% linspace can return 8.000000000000001e-07 for a nominal 800 ns, so a bare
% 'tau > cpDur' misclassifies it. Guard with a relative tolerance.
powerBeyondCP = sum(pdp(tau > cpDur*(1+1e-9)));   % fraction of power -> ISI
if     powerBeyondCP == 0,     isiSeverity = 'none';
elseif powerBeyondCP < 0.05,   isiSeverity = 'mild';
else,                          isiSeverity = 'severe';
end

% ---------------- correlation square roots ----------------
RtxH = sqrtm(corrMatrix(Ntx,rhoTx));
RrxH = sqrtm(corrMatrix(Nrx,rhoRx));

% ---------------- shape / symbol count ----------------
switch lower(cfg.domain)
    case 'time'
        if isvector(txSignal), txSignal = txSignal(:); end
        nSamp = size(txSignal,1);
        nSym  = floor(nSamp/(cfg.nFFT+cfg.nCP));
        assert(nSym >= 1, 'txSignal shorter than one OFDM symbol.');
        nGain = nSamp;  gainRate = sr;
    case 'grid'
        if ismatrix(txSignal)
            txSignal = reshape(txSignal,size(txSignal,1),size(txSignal,2),1);
        end
        nSym  = size(txSignal,2);
        nGain = nSym;   gainRate = 1/Tsym;
    otherwise
        error('cfg.domain must be ''time'' or ''grid''.');
end
% holdOverSymbols freezes the channel: one gain sample reused everywhere.
if cfg.holdOverSymbols, nGain = 1; end

% ============================================================================
% TIME-VARYING TAP GAINS. Each tap has its OWN Jakes process, so Doppler is
% genuine fading rather than one frequency offset shared by every tap.
% Validated: unit power, autocorrelation matches J0(2*pi*fd*tau).
% ============================================================================
% The LOS ray is DETERMINISTIC but NOT static: a moving terminal shifts it by
% a Doppler frequency. Freezing it (los = const) is only correct for a
% stationary link. For a fast terminal it removes the LOS contribution to ICI
% entirely - which for the LEO satellite (K=1) meant HALF the channel power
% produced no ICI, understating the error floor by roughly 2x.
% A strong LOS therefore makes ICI WORSE for a fast terminal, not better.
losPhase = exp(1i*2*pi*fd*((0:nGain-1).'/gainRate));   % [nGain x 1]

g = zeros(Nrx,Ntx,nTaps,nGain);
for t = 1:nTaps
    Graw = zeros(Nrx,Ntx,nGain);
    for rxi = 1:Nrx
        for txi = 1:Ntx
            Graw(rxi,txi,:) = jakesTap(nGain, fd, gainRate);
        end
    end
    if cfg.K > 0 && t == 1                       % rank-1 LOS on the first tap
        for n = 1:nGain
            Graw(:,:,n) = sqrt(cfg.K/(cfg.K+1))*losPhase(n)*ones(Nrx,Ntx) ...
                        + sqrt(1/(cfg.K+1))*Graw(:,:,n);
        end
    end
    for n = 1:nGain                              % correlation + PDP scaling
        g(:,:,t,n) = sqrt(pdp(t)) * (RrxH * Graw(:,:,n) * RtxH);
    end
end

% ============================================================================
% NOISE VARIANCE - from the MEASURED TRANSMITTED power (see header).
% ============================================================================
txPow    = mean(abs(txSignal(:)).^2);
sigPow   = txPow * Ntx;                 % mean received power per Rx antenna
noiseVar = sigPow / snrToLinear(cfg);

% ============================================================================
% APPLY THE CHANNEL
% ============================================================================
switch lower(cfg.domain)

  case 'time'
    assert(size(txSignal,2)==Ntx,'txSignal columns must equal Ntx.');
    rxFaded = zeros(nSamp,Nrx);
    for rxi = 1:Nrx
        acc = zeros(nSamp,1);
        for txi = 1:Ntx
            for t = 1:nTaps
                d = delaySamp(t);
                if d >= nSamp, continue; end
                gt = gainVec(g,rxi,txi,t,nGain,nSamp);     % [nSamp x 1]
                acc(1+d:nSamp) = acc(1+d:nSamp) ...
                    + gt(1+d:nSamp).*txSignal(1:nSamp-d,txi);
            end
        end
        rxFaded(:,rxi) = acc;
    end
    rxSignal = rxFaded;
    if cfg.addNoise
        rxSignal = rxSignal + sqrt(noiseVar/2)*(randn(nSamp,Nrx)+1i*randn(nSamp,Nrx));
    end

    % --- H: per-symbol AVERAGE channel (Stage 4 convention). Within-symbol
    %     variation is deliberately NOT captured, so it leaks through as ICI.
    %     Vectorised: gAvg [nTaps x nSym], then H = E * gAvg with the DFT
    %     kernel E [nFFT x nTaps]. Fast enough for BER sweeps.
    symLen = cfg.nFFT + cfg.nCP;
    E = exp(-1i*2*pi*(0:cfg.nFFT-1).'*delaySamp(:).'/cfg.nFFT);   % [nFFT x nTaps]
    H4 = zeros(Nrx,Ntx,cfg.nFFT,nSym);
    for rxi = 1:Nrx
        for txi = 1:Ntx
            gAvg = zeros(nTaps,nSym);
            for t = 1:nTaps
                gt = gainVec(g,rxi,txi,t,nGain,nSamp);
                blk = reshape(gt(1:nSym*symLen), symLen, nSym);
                gAvg(t,:) = mean(blk,1);          % average over each symbol
            end
            H4(rxi,txi,:,:) = E * gAvg;           % [nFFT x nSym]
        end
    end

  case 'grid'
    [nSc,~,txAnt] = size(txSignal);
    assert(txAnt==Ntx,'txSignal 3rd dim must equal Ntx.');
    assert(nSc <= cfg.nFFT,'nSc cannot exceed nFFT.');
    cfg = setdef(cfg,'scIdx',1:nSc);
    assert(numel(cfg.scIdx)==nSc,'numel(cfg.scIdx) must equal nSc.');

    rxFaded = zeros(nSc,nSym,Nrx);
    H4      = zeros(Nrx,Ntx,nSc,nSym);
    for s = 1:nSym
        gi = min(s,nGain);                      % frozen if holdOverSymbols
        for ii = 1:nSc
            f = cfg.scIdx(ii);
            Hsc = zeros(Nrx,Ntx);
            for t = 1:nTaps
                Hsc = Hsc + g(:,:,t,gi)*exp(-1i*2*pi*(f-1)*delaySamp(t)/cfg.nFFT);
            end
            H4(:,:,ii,s) = Hsc;
            x = reshape(txSignal(ii,s,:),Ntx,1);
            rxFaded(ii,s,:) = reshape(Hsc*x,1,1,Nrx);
        end
    end
    rxSignal = rxFaded;
    if cfg.addNoise
        rxSignal = rxSignal + sqrt(noiseVar/2)* ...
            (randn(nSc,nSym,Nrx)+1i*randn(nSc,nSym,Nrx));
    end
end

% ---------------- pack H (drop symbol dim when frozen) ----------------
if cfg.holdOverSymbols
    H3 = H4(:,:,:,1);                                  % [Nrx x Ntx x nSc/nFFT]
    if isMimo, H = H3; else, H = reshape(H3(1,1,:),[],1); end
else
    if isMimo, H = H4;                                 % [Nrx x Ntx x nSC x nSym]
    else,      H = reshape(H4(1,1,:,:), size(H4,3), nSym);
    end
end

% ---------------- dual-limit metrics ----------------
Bc       = 1/(5*max(rms_s,eps));
scSpace  = sr/cfg.nFFT;
BcOverSC = Bc/scSpace;

% J0 0.5-crossing, matching Stage 4. J0(x)=0.5 at x = 1.5211.
if fd < 1e-6, Tc = Inf; else, Tc = 1.5211/(2*pi*fd); end
TcOverTsym = Tc/Tsym;

if strcmp(isiSeverity,'severe'),   bind = 'FREQUENCY (severe ISI)';
elseif TcOverTsym < 1,             bind = 'TIME (ICI, Doppler)';
elseif strcmp(isiSeverity,'mild'), bind = 'FREQUENCY (mild ISI)';
else,                              bind = 'none (comfortable)';
end

info = struct( ...
    'noiseVariance',noiseVar,'measuredTxPow',txPow,'meanRxPow',sigPow, ...
    'delayTaps_s',tau,'pdp',pdp, ...
    'rmsDelaySpread_s',rms_s,'maxDelay_s',maxDelay,'cpDuration_s',cpDur, ...
    'cpExceeded',cpExceeded,'powerBeyondCP',powerBeyondCP, ...
    'isiSeverity',isiSeverity,'bindingLimit',bind, ...
    'coherenceBW_Hz',Bc,'BcOverSC',BcOverSC, ...
    'fd_Hz',fd,'coherenceTime_s',Tc,'TcOverTsym',TcOverTsym, ...
    'symbolDuration_s',Tsym,'K',cfg.K,'speed_kmh',cfg.speed_kmh, ...
    'rhoTx',rhoTx,'rhoRx',rhoRx,'domain',cfg.domain,'snrMode',cfg.snrMode, ...
    'holdOverSymbols',cfg.holdOverSymbols,'numSymbols',nSym);

% ---------------- diagnostics ----------------
if cfg.verbose
    fprintf('---------------- ofdm_channel ----------------\n');
    if isfield(cfg,'environment'), fprintf('  environment    : %s\n',cfg.environment); end
    fprintf('  domain         : %s   | antennas %dx%d (Tx x Rx)\n',cfg.domain,Ntx,Nrx);
    if isMimo, fprintf('  correlation    : rhoTx=%.2f rhoRx=%.2f\n',rhoTx,rhoRx); end
    fprintf('  RMS delay      : %.0f ns | max delay %.0f ns | CP %.0f ns\n', ...
        rms_s*1e9, maxDelay*1e9, cpDur*1e9);
    fprintf('  CP exceeded    : %s | power beyond CP = %.2f%% -> %s ISI\n', ...
        ternary(cpExceeded,'YES','no'), powerBeyondCP*100, isiSeverity);
    fprintf('  FREQUENCY      : Bc = %.0f kHz | Bc/SC = %.2f\n',Bc/1e3,BcOverSC);
    fprintf('  TIME           : fd = %.0f Hz | fd*Tsym = %.4f | Tc/Tsym = %.2e\n', ...
        fd, fd*Tsym, TcOverTsym);
    fprintf('  BINDING LIMIT  : %s\n',bind);
    fprintf('  SNR mode       : %s (%.1f dB)\n',cfg.snrMode,cfg.snrdB);
    fprintf('  measured txPow : %.6f  (MEASURED, not assumed = 1)\n',txPow);
    fprintf('  noise variance : %.4e\n',noiseVar);
    fprintf('  holdOverSymbols: %s (%d symbols)\n', ...
        ternary(cfg.holdOverSymbols,'true (frozen)','false (Doppler varies)'), nSym);
    fprintf('----------------------------------------------\n');
end
end


% ============================================================================
% expand a (possibly frozen) tap gain to a full-length time vector
% ============================================================================
function gt = gainVec(g,rxi,txi,t,nGain,nSamp)
    if nGain == 1
        gt = g(rxi,txi,t,1)*ones(nSamp,1);
    else
        gt = reshape(g(rxi,txi,t,:),[],1);
    end
end

% ============================================================================
function snrLin = snrToLinear(cfg)
    switch lower(cfg.snrMode)
        case 'snr'
            snrdB = cfg.snrdB;
        case 'ebno'
            snrdB = cfg.snrdB + 10*log10(cfg.bitsPerSym) ...
                              + 10*log10(cfg.nDataSC/cfg.nFFT);
        otherwise
            error('cfg.snrMode must be ''snr'' or ''ebno''.');
    end
    snrLin = 10^(snrdB/10);
end

% ============================================================================
% Sum-of-sinusoids Jakes tap gain (unit power, time-varying).
% fd = 0 gives a static Rayleigh gain held over the frame (block fading).
% Validated: mean|g|^2 = 1.00; R(tau) matches J0(2*pi*fd*tau) to <0.002.
% ============================================================================
function gt = jakesTap(n, fd, rate)
    if fd <= 0 || n == 1
        gt = ((randn + 1i*randn)/sqrt(2))*ones(n,1);
        return;
    end
    M = 16;
    t = (0:n-1).'/rate;
    m = 1:M;
    theta = 2*pi*rand - pi;
    alpha = (2*pi*m - pi + theta)/(4*M);
    phi = 2*pi*rand(1,M) - pi;
    psi = 2*pi*rand(1,M) - pi;
    arg = 2*pi*fd*(t*cos(alpha));               % [n x M]
    gi = sqrt(2/M)*sum(cos(arg + phi),2);
    gq = sqrt(2/M)*sum(cos(arg + psi),2);
    gt = (gi + 1i*gq)/sqrt(2);
end

% ============================================================================
function p = envPreset(name)
    switch lower(name)
        % RMS delay spreads anchored to published channel models:
        %   indoor  -> IEEE 802.11n TGn (B=15, D=50, F=150 ns)
        %   urban   -> dense urban, COST207 Typical-Urban range
        %   highway -> ITU-R M.1225 Vehicular-A (370 ns)
        %   sat     -> 3GPP TR 38.811 NTN (delay ~100 ns; 130 kHz Doppler @5.2GHz,
        %              = the S-band 48 kHz@2GHz worst case for a 7.5 km/s LEO)
        case 'indoor_home',   p = struct('rms',15e-9,  'K',8.0,'v',3);   % TGn-B residential
        case 'indoor_office', p = struct('rms',50e-9,  'K',3.0,'v',3);   % TGn-D office
        case 'hotspot',       p = struct('rms',150e-9, 'K',1.0,'v',5);   % TGn-F large space
        case 'urban',         p = struct('rms',700e-9, 'K',0.5,'v',50);  % dense urban
        case 'highway',       p = struct('rms',370e-9, 'K',0.3,'v',120); % ITU Veh-A, open V2V
        case 'satellite',     p = struct('rms',100e-9, 'K',1.0,'v',27000); % 3GPP NTN S-band
        otherwise
            error(['ofdm_channel: unknown environment "%s". Use indoor_home | ' ...
                   'indoor_office | hotspot | urban | highway | satellite'],name);
    end
end

function R = corrMatrix(n,rho)
    idx = (1:n).';
    R = rho.^abs(idx-idx.');
    R = (R+R')/2;
end

function s = setdef(s,f,v)
    if ~isfield(s,f) || isempty(s.(f)), s.(f) = v; end
end

function out = ternary(cond,a,b)
    if cond, out = a; else, out = b; end
end
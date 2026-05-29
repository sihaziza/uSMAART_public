function [signalE_aligned, eTTL_aligned, tO, tE_corrected, model] = alignClockDriftFromTTL( ...
    eTTL, oTTL, signalE, FsE, FsO, method, edgeType)
% ALIGNCLOCKDRIFTFROMTTL
%
% Corrects clock drift between two acquisition systems using TTL time series.
%
% INPUTS
%   eTTL    : TTL square-wave time series from E system
%   oTTL    : TTL square-wave time series from O system
%   signalE : signal recorded on E system
%   FsE     : sampling rate of E system
%   FsO     : sampling rate of O system
%   method  : 'polyfit' or 'fitlm'
%   edgeType: 'rising', 'falling', or 'both'
%
% OUTPUTS
%   signalE_aligned : signalE corrected/resampled onto O time base
%   tO              : O-system time vector
%   model           : alignment model and diagnostics

if nargin < 6 || isempty(method)
    method = 'fitlm';
end

if nargin < 7 || isempty(edgeType)
    edgeType = 'rising';
end

%% Detect TTL event times
ttlE_times = find(diff(eTTL)>0)/FsE;
ttlO_times = find(diff(oTTL)>0)/FsO;

% ttlE_times = detectTTLEdges(oTTL, FsO, edgeType);
% ttlO_times = detectTTLEdges(oTTL, FsO, edgeType);

%% Match number of TTL events
nTTL = min(numel(ttlE_times), numel(ttlO_times));

ttlE_times = ttlE_times(1:nTTL);
ttlO_times = ttlO_times(1:nTTL);

if nTTL < 3
    error('Not enough matched TTL events to estimate clock drift.');
end

%% Fit affine time transform
% Goal: tO ≈ a*tE + b

switch lower(method)

    case 'polyfit'
        p = polyfit(ttlE_times, ttlO_times, 1);
        a = p(1);
        b = p(2);

        model.type = 'polyfit';

    case 'fitlm'
        mdl = fitlm(ttlE_times(:), ttlO_times(:), ...
            'linear', 'RobustOpts', 'on');

        b = mdl.Coefficients.Estimate(1);
        a = mdl.Coefficients.Estimate(2);

        model.type = 'fitlm';
        model.mdl = mdl;

    otherwise
        error('Unknown method. Use ''polyfit'' or ''fitlm''.');
end

%% Build corrected E timebase
nE = size(signalE,1);
tE = (0:nE-1)' / FsE;

% tE_corrected = a * tE; % not including the offset + b;
tE_corrected = a * tE + b;
% sum(isnan(tE_corrected))

%% Build O timebase
nO = size(oTTL,1);
tO = (0:nO-1)' / FsO;

%% Resample signalE onto O timebase

% for analogue channel
for iChannel=1:size(signalE,2)
signalE_aligned(:,iChannel) = interp1( ...
    tE_corrected, ...
    signalE(:,iChannel), ...
    tO, ...
    'pchip');
end

% sum(isnan(signalE_aligned))

% for TTL channel
eTTL_aligned = interp1(tE_corrected,eTTL,tO,'next','extrap');
% sum(isnan(eTTL_aligned))

%% Diagnostics
ttlE_predicted_in_O_time = a * ttlE_times + b;
residualLag = ttlO_times - ttlE_predicted_in_O_time;

model.a = a;
model.b = b;
model.ppm = (a - 1) * 1e6;
model.ttlE_times = ttlE_times;
model.ttlO_times = ttlO_times;
model.residualLag = residualLag;
model.residualLag_ms = residualLag * 1000;
model.edgeType = edgeType;
model.nTTL = nTTL;

fprintf('\n=== CLOCK DRIFT ALIGNMENT ===\n');
fprintf('Method            : %s\n', method);
fprintf('TTL edge type     : %s\n', edgeType);
fprintf('Matched TTLs      : %d\n', nTTL);
fprintf('Scale factor a    : %.12f\n', a);
fprintf('Offset b          : %.6f s\n', b);
fprintf('Clock drift       : %.3f ppm\n', model.ppm);
fprintf('Residual lag RMS  : %.3f ms\n', rms(model.residualLag_ms));

end

%% ------------------------------------------------------------------------
function ttlTimes = detectTTLEdges(ttl, Fs, edgeType)

ttl = ttl(:);

% Robust threshold for irregular TTL amplitude
lo = prctile(ttl, 5);
hi = prctile(ttl, 95);
thr = (lo + hi) / 2;

ttlBinary = ttl > thr;

dTTL = diff(ttlBinary);

switch lower(edgeType)

    case 'rising'
        edgeIdx = find(dTTL == 1) + 1;

    case 'falling'
        edgeIdx = find(dTTL == -1) + 1;

    case 'both'
        edgeIdx = find(dTTL ~= 0) + 1;

    otherwise
        error('edgeType must be ''rising'', ''falling'', or ''both''.');
end

ttlTimes = (edgeIdx - 1) / Fs;

end
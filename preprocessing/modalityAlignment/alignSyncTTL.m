function [range_ttlE,range_ttlO]=alignSyncTTL(ttlE,ttlO,fps,varargin)
% align two time trace based off a common syncTTL. Assume that both are at
% the same sampling rate
% [ttlE_cal,ttlO_cal]=alignRecording(ttlE,ttlO,fps)

% 20211209 SH - fixing bug when length after shift does not match anymore

%% OPTIONS
options.savePath=[];
options.plotFigure=true;

%% UPDATE OPTIONS
if nargin>3
    options=getOptions(options,varargin);
end

%%
if ~iscolumn(ttlE)
    ttlE=ttlE';
end

if ~iscolumn(ttlO)
    ttlO=ttlO';
end

% normalize [0 1] for analogue traces
ttlE=rescale(ttlE,0,1);ttlE=ttlE-median(ttlE);
ttlO=rescale(ttlO,0,1);ttlO=ttlO-median(ttlO);

% Find the time lag between using xcorr
[r,lags]=xcorr(ttlE,ttlO);
[~,idShift]=max(r);
shift=lags(idShift);

if shift>0
    disp('ePhys is in advance')
    
    range_ttlO=1:length(ttlO);
    range_ttlE=shift+1:min(length(ttlE),length(ttlO)+shift);
    range_ttlO=range_ttlO(1:length(range_ttlE));
    
    ttlO_cal=ttlO(range_ttlO);
    ttlE_cal=ttlE(range_ttlE);
    
else
    disp('TEMPO is in advance')
    shift=abs(shift);
    range_ttlO=shift+1:length(ttlE)+shift;
    range_ttlE=1:length(ttlE);
    
    ttlO_cal=ttlO(range_ttlO);
    ttlE_cal=ttlE(range_ttlE);
end

if options.plotFigure
    
    figHandle=figure('Name','SyncTTL alignment');
    lnWidth=1.5;
    
    subplot(2,3,1)
    plot(lags/fps,r)
    hold on
    plot([0 0],[min(r) max(r)],'k')
    xlim([-30 30])
    title('Time lag - pre alignment')
    
    subplot(2,3,[2 3])
    plot(getTime(ttlE,fps),ttlE-1,'linewidth',lnWidth)
    hold on
    plot(getTime(ttlO,fps),ttlO,'linewidth',lnWidth)
    hold off
    title('Raw data in their own time frame')
    
    subplot(2,3,4)
    [r,lags]=xcorr(ttlE_cal,ttlO_cal);
    
    plot(lags/fps,r)
    hold on
    plot([0 0],[min(r) max(r)],'k')
    xlim([-30 30])
    xlabel('Time (s)')
    title('Time lag - post alignment')
    
    time=getTime(ttlO_cal,fps);
    
    subplot(2,3,[5 6])
    plot(time,ttlE_cal-1,'linewidth',lnWidth)
    hold on
    plot(time,ttlO_cal,'linewidth',lnWidth)
    hold off
    xlabel('Time (s)')
    title('Multimodal re-alignment')
    
    if options.savePath
        savePDF(figHandle,'alignSyncTTL',options.savePath)
    end
end
end
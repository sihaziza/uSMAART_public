function [output]=rasterERP(wave,ttl,fs,varargin)
% e.g. rasterERP(lfp,ttl,fs)
% options.cRange='auto'; % usually standard deviation
% options.title='ERP raster';
% options.stimBand=[];
% options.figHandle=[];
% options.figure=true;
% options.ttlIDX=[];
% options.getShuffle=false;
% options.preNormalize='median'; % mean %median % none
% options.stimLength=[];
% options.baselinePrePost=1; % 1sec baseline

%% DEFAULT OPTIONS
options.cRange='auto'; % usually standard deviation
options.title='ERP raster';
options.stimBand=[];
options.figHandle=[];
options.figure=true;
options.ttlIDX=[];
options.getShuffle=false;
options.preNormalize='mean'; % mean %median % none %zcore
options.stimLength=[];
options.baselinePrePost=1; % 1sec baseline

%% USER-DEFINED INPUT OPTIONS
if nargin>3
    options=getOptions(options,varargin);
end

%% Plot event-related power change

[output]=getStimEpoch(wave,ttl,fs,...
    'stimLength',options.stimLength,'baselinePrePost',options.baselinePrePost,'preNormalize',options.preNormalize,'getShuffle',options.getShuffle);

array=output.arrayRaw;
stimBand=output.stimBand;

%%
n=size(array,2);
t=getTime(array,fs)+stimBand(1);
output.time=getTime(array,fs)+stimBand(1);

if options.figure
    
if ~isempty(options.figHandle)
    figure(options.figHandle)
else
    figure('DefaultAxesFontSize',14,'color','w');
end
    
imagesc(t,0:n,array')
shading 'flat' %'interp'
caxis(options.cRange)
colormap('redblue')
colorbar('Location','southoutside');
% c.Label.String = 'SD norm.';
xlim([t(1) t(end)])
hold on
plot([0 0],[0 n],'--k','LineWidth',1)
hold off
title(options.title)
ylabel('Event #')
axis tight

% figure('DefaultAxesFontSize',14,'color','w');
% [L,H]=bounds(mean(array,2),'all');
% yRange=[L H]*1.1;
% plotErrorBar1(array,'x_axis',t,'error','sem');
% ylabel('{\Delta}F/F (%)')
% hold on
% plot([0 0],yRange,'--k','LineWidth',1)
% hold off
% xlabel('Time (s)')
% ylabel('SD Norm.')
% axis tight
end

end
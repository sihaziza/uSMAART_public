function [outputParse,f]=powerERP(wave,ttl,fs,varargin)
% EXAMPLE: [[outputParse,f]=powerERP(wave,ttl,fs);
% options.figHandle=[];
% options.figure=true;
% options.ttlIDX=[];
% options.baselinePrePost=2; % 1sec baseline
% options.getShuffle=false;
% options.preNormalize='median'; % mean %median % none
% options.stimLength=1;
% options.BW=[1 150];

% DEFAULT Options
options.figHandle=[];
options.figure=true;
options.ttlIDX=[];
options.baselinePrePost=2; % 1sec baseline
options.getShuffle=false;
options.preNormalize='none'; % mean %median % none
options.stimLength=1;
options.BW=[1 150];

%% USER-DEFINED INPUT OPTIONS
if nargin>3
    options=getOptions(options,varargin);
end

%%

BW=options.BW;

% compute sensory-stim evoked response on temporal and frequency domaine
[wt,f] = cwt(wave,fs,'FrequencyLimits',[0.5 floor(fs/2)]);
disp('continuous wavelet transform done...')

%% Plot event-related power change

[outputParse]=getStimEpoch(abs(wt)',ttl,fs,...
    'stimLength',options.stimLength,'baselinePrePost',options.baselinePrePost,'preNormalize',options.preNormalize,'getShuffle',false);
outputParse.f=f;

if options.figure
    if ~isempty(options.figHandle)
        figure(options.figHandle)
    else
        fig=figure('DefaultAxesFontSize',14,'color','w');
    end
    
    time=getTime(outputParse.arrayRaw,fs)+outputParse.stimBand(1);
    
    if isempty(options.ttlIDX)
        avgPow=squeeze(mean(outputParse.arrayRaw,2));
    else
        disp('keeping only a subset of events, based off user input')
        avgPow=squeeze(mean(outputParse.arrayRaw(:,options.ttlIDX,:),2));
    end
    
%     %%
% tmp1 = abs(avgPow');
% t1 = size(tmp1,2);
% tmp1 = tmp1';
% minv = min(tmp1,[],1);
% tmp1 = (tmp1-minv(ones(1,t1),:));
% 
% maxv = max(tmp1);
% maxvArray = maxv(ones(1,t1),:);
% indx = maxvArray<eps;
% tmp1 = 100*(tmp1./maxvArray);
% tmp2 = 1+fix(tmp1);
% tmp2(indx) = 1;
% tmp2 = tmp2';
% 
% pcolor(time,f,tmp2);shading interp
% colormap(jet)
% colorbar
% 
%%
% epoch=outputParse.stimBand;
% 
% tmp1 = abs(avgPow);
% baseline = median(tmp1(1:round(1*fs),:),1);
% tmp1 = 100*((tmp1-baseline)./baseline)';
% 
% pcolor(time,f,tmp1);shading interp
% colormap(jet)
% colorbar
% ylim([1 55])
% caxis([0 300])
% imagesc(time,f,tmp1);

    %%
    imagesc(time,f,(avgPow)');%,'CDataMapping','scaled');
    ax=gca();
    % chigh=round(max(avgPow,[],'all'),3);
    % chigh=2;
    ax.YLim = BW;
    ax.XLim = [min(time) max(time)];
    % ax.CLim = [0 chigh];
    ax.Layer = 'top';
    ax.YDir = 'normal';
    ax.YScale = 'log';
    ax.YTick=[1 2 5 10 20 50 100];
    ax.YTickLabel={'1' '2' '5' '10' '20' '50' '100'};
    xlabel('Time (s)')
    ylabel('Frequency (Hz)')
    c=colorbar;
    c.Label.String='Signal magnitude (%)';
    title('power ERP')
    hold on
    plot([0 0],ax.YLim,'--k')
    hold off
end
end
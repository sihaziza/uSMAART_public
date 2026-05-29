function [allSpindles,ttlSpindles]=findSpindles(LFP,varargin)

%% DEFAULT Options
options.verbose=true;
options.plotFigure=false;
options.figHandle=[];
options.plotRaster=false;
options.saveFig=false;
options.plotAverage=false;
options.CSD=false;

options.samplingRate=2000;
options.BW=[15 20];
options.window=250; % in ms

options.minPeakProm=20; % in percent
options.minPeakDist=1000; % in ms
options.minPeakWidth=100; % in ms
options.maxPeakWidth=2000; % in ms

%%
% USER-DEFINED INPUT OPTIONS
if nargin>1
    options=getOptions(options,varargin);
end

%% STEP 1: coarse - identify likely ripple localization

fs=options.samplingRate;
BW=options.BW;
% 
% [sw,f] =cwt(LFP,fs,'FrequencyLimits',[1 100],'amor');
% 
% id_on=find(f>BW(1),1,'first');
% id_off=find(f<BW(2),1,'last');
%   
% spinPW_Trace=rescale(mean(abs(sw(id_on:id_off,:))),0,1);

spindles=bpFilter1D(LFP,BW,fs);
rippleBand=spindles.^2;
rippleBand=movmax(rippleBand,fs/50);
rippleBand=smooth(rippleBand,fs/50);
rippleBand=rescale(rippleBand,0,1);

[~,locs]=findpeaks(rippleBand,...
    'MinPeakProminence',options.minPeakProm/100,...
    'MinPeakDistance',options.minPeakDist/1000*fs,...
    'MinPeakWidth',options.minPeakWidth/1000*fs,...
    'MaxPeakWidth',options.maxPeakWidth/1000*fs);

allSpindles=[];tempRange=[];range=[];
nSpin=numel(locs);
ttlSpindles=zeros(size(spindles));
periT=250; % in ms
window=[-250 250]; % in ms
for iSpin=1:nSpin
    
    tempRange=locs(iSpin)-fs*periT/1000:locs(iSpin)+fs*periT/1000-1;
    
    % find min value to align all ripples to each other
    temp=rescale((spindles(tempRange)),-1,1); 
    %     [~,id]=  findpeaks(-temp,'NPeaks', 1,'MinPeakHeight',0.5,'MinPeakDistance',10,'MaxPeakWidth',8,'MinPeakWidth',5);
    [~,id]=min(temp);
    
    ttlSpindles(locs(iSpin)-fs*periT/1000+id)=1;
    
    range(:,iSpin)=locs(iSpin)+(id-fs*periT/1000)+fs*window(1)/1000:locs(iSpin)+(id-fs*periT/1000)+fs*window(2)/1000-1;
    
    allSpindles(:,iSpin)=LFP(range(:,iSpin));
end

%%

wave=bpFilter1D(allSpindles,[5 inf],fs);
time=getTime(wave,fs)*1000-periT;

if options.plotRaster
    figure
    imagesc(time,1:nSpin,wave')
end

if options.plotAverage
    plotErrorBar1(wave,'x_axis',time);
end

disp(['found ' num2str(nSpin) ' spindles. Let"s clean this up'])

end
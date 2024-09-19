function [cleanedRipple,ttlSWR]=findRipples(LFP,varargin)
% heuristic detection of ripples
% example: [allRipples,range,ttlSWR]=findRipples(LFP);
% example: [allRipples,range,ttlSWR]=findRipples(LFP,'samplingRate',1000,'BW',[120 250]);

%% DEFAULT Options
options.verbose=true;
options.plotFigure=false;
options.figHandle=[];
options.plotRaster=false;
options.saveFig=false;
options.plotAverage=false;
options.CSD=false;

options.samplingRate=2000;
options.BW=[120 200];
options.window=100; % in ms

options.minPeakProm=10; % in percent
options.minPeakDist=10; % in ms
options.minPeakWidth=10; % in ms
options.maxPeakWidth=200; % in ms

%%
% USER-DEFINED INPUT OPTIONS
if nargin>1
    options=getOptions(options,varargin);
end

%% STEP 1: coarse - identify likely ripple localization

fs=options.samplingRate;
BW=options.BW;

ripple=bpFilter1D(LFP,BW,fs);
rippleBand=ripple.^2;
rippleBand=movmax(rippleBand,fs/50);
rippleBand=smooth(rippleBand,fs/50);
rippleBand=rescale(rippleBand,0,1);

[~,locs]=findpeaks(rippleBand,...
    'MinPeakProminence',options.minPeakProm/100,...
    'MinPeakDistance',options.minPeakDist/1000*fs,...
    'MinPeakWidth',options.minPeakWidth/1000*fs,...
    'MaxPeakWidth',options.maxPeakWidth/1000*fs);

allRipples=[];
nSWR=numel(locs);
ttlSWR=zeros(size(ripple));
periT=25; % in ms
window=[-150 150]; % in ms
for iSWR=1:nSWR
    
    tempRange=locs(iSWR)-fs*periT/1000:locs(iSWR)+fs*periT/1000-1;
    
    % find min value to align all ripples to each other
    temp=rescale((ripple(tempRange)),-1,1);
    %     [~,id]=  findpeaks(-temp,'NPeaks', 1,'MinPeakHeight',0.5,'MinPeakDistance',10,'MaxPeakWidth',8,'MinPeakWidth',5);
    [~,id]=min(temp);
    
    ttlSWR(locs(iSWR)-fs*periT/1000+id)=1;
    
    range(:,iSWR)=locs(iSWR)+(id-fs*periT/1000)+fs*window(1)/1000:locs(iSWR)+(id-fs*periT/1000)+fs*window(2)/1000-1;
    
    allRipples(:,iSWR)=LFP(range(:,iSWR));
    %     allRipples32(:,:,iSWR)=LFP(range(:,iSWR),:);
end

%% STEP 2 : refine - clean up bad ripples and identify start-end of ripple

% figure position for home workstation [1383,953,1920,970] / lab workstation [6,1241,1920,970]
% surfacePro [77,148,1262,676]
options.positionFig=[77,148,1262,676];

figure('defaultaxesfontsize',12,'color','w','Position',options.positionFig);
k=1;
for iSWR=1:size(allRipples,2)
    wave=bpFilter1D(allRipples(:,iSWR),[5 inf],fs);
    
    [sw,f] =cwt(wave',fs,'FrequencyLimits',options.BW,'amor');
    % ylim([0 0.30])
    
    id_on=find(f>135,1,'first');
    id_off=find(f<200,1,'last');
    
    temp=rescale(mean(abs(sw(id_on:id_off,:))),0,1);
    temp=temp(151:end-150);
    
    % plot(temp)
    options.minPeakProm=30; % in percent
    options.minPeakDist=20; % in ms
    options.minPeakWidth=10; % in ms
    options.maxPeakWidth=200; % in ms
    options.threshold=30; % in percent
    
    [pks,locs,~]=findpeaks(temp,...
        'MinPeakProminence',options.minPeakProm/100,...
        'MinPeakDistance',options.minPeakDist/1000*fs,...
        'MinPeakWidth',options.minPeakWidth/1000*fs,...
        'MaxPeakWidth',options.maxPeakWidth/1000*fs,...
        'Annotate','extents','WidthReference','halfprom');
    
    thres=options.threshold/100; % 50% of the max = FWHM
    
    time=getTime(wave,fs);
    waveEpoch=[];timeEpoch=[];
    rippleStartTime=[];duration=[];
    for iLoc=1%:numel(locs)
        
        idx_on=max(locs(iLoc)-0.1*fs,1);
        idx_off=min(locs(iLoc)+0.1*fs,300);
        
        % start from the center
        idx_on=fliplr(temp(idx_on:locs(iLoc)));plot(idx_on)
        idx_on=find(idx_on<thres*pks(iLoc),1,'first');
        if isempty(idx_on)
            idx_on=1;
        end
        idx_on=max(locs(iLoc)-idx_on,1);
        
        idx_off=temp(locs(iLoc):idx_off);
        idx_off=find(idx_off<thres*pks(iLoc),1,'first');
        if isempty(idx_off)
            idx_off=length(wave);
        end
        idx_off=min(locs(iLoc)+idx_off,length(wave));
        
        %     rippleStartTime(iLoc)=idx_on;
        %     duration(iLoc)=idx_off-idx_on; % in ms
        
        waveEpoch=[waveEpoch; wave(idx_on+150:idx_off+149)];
        timeEpoch=[timeEpoch time(idx_on+150:idx_off+149)];
    end
    
    
    plot(time,wave,'k')
    hold on
    plot(timeEpoch,waveEpoch,'r')
    plot(timeEpoch(1),waveEpoch(1),'*g')
    plot(timeEpoch(end),waveEpoch(end),'*g')
    % find first min peak to align all ripples to each other
    temp=bpFilter1D(waveEpoch,[100 inf],fs);
    [~,id_pk]=min(temp(1:round(2*fs/200))); %find min over 1 oscillation from start
    plot(timeEpoch(id_pk),waveEpoch(id_pk),'*g')
    hold off
    
    opts.Interpreter = 'tex';
    % Include the desired Default answer
    opts.Default = 'No';
    % Use the TeX interpreter to format the question
    quest = 'keep this ripple?';
    answer = questdlg(quest,'ripple clean-up','Yes','No',opts);
    
    if strcmpi(answer,'yes')
        cleanedRipple(k).start=timeEpoch(1);
        cleanedRipple(k).minPks=timeEpoch(id_pk);
        cleanedRipple(k).end=timeEpoch(end);
        cleanedRipple(k).index=iSWR;
        k=k+1;
    end
    
end
%%
time=getTime(allRipples,fs)*1000-100;

if options.plotRaster
    figure
    imagesc(time,1:nSWR,allRipples')
end

if options.plotAverage
    plotErrorBar1(allRipples,'x_axis',time);
end

disp(['found ' num2str(nSWR) ' ripples. Let"s clean this up'])

if options.CSD
    CSDoutput=[];
    for i=1:65
        [CSDoutput(:,:,i)]  = compCSD(rippleLFP(:,:,i)/1000,fs,50*1e-6);
    end
    
    imagesc(mean(CSDoutput,3)')
end
end
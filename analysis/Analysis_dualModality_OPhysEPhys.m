%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% ANALYSIS of mPFC PV-INs uSMAART/LFP data %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
% close all

% pick mouse number 1 or 2
mouseN=2;

% !!!change this line!!!
% path=['yourPath to uSMAART .mat & LFP .rhd data']; % !!!change this line!!!

path=['Z:\Mie\m' num2str(mouseN)];

file=dir(fullfile(path,'**/*.mat'));
allData=[];

nbytes = fprintf('processing %d%% \n', 0);
nFile=numel(file);
for iFile=1:nFile
    [temp_pre,meta]=loaduSMAART2mat(fullfile(file(iFile).folder,file(iFile).name),'figure',false,'verbose',false);
    allData=[allData; temp_pre];
    
    fprintf(repmat('\b',1,nbytes))
    nbytes = fprintf('processing %d%% \n', round(iFile/nFile*100));
end

[u,~]=find(isnan(allData));
allData=allData(setdiff(1:length(allData),unique(u)),:);

fs=meta.fs;

% T=getTime(allData,fs);
% figure; plot(T(1:100:end)./60,allData(1:100:end,1:5))

oTTL=allData(:,6:7);
signal=allData(:,3:5);

T=getTime(signal,fs);
figure; plot(T(1:2:end)./60,signal(1:2:end,:))

%% plot photobleaching time course

figure('DefaultAxesFontSize',14,'color','w');
for i=1:3
    temp=signal(:,i);temp=(temp)./max(temp(1:60*fs));
    T=getTime(temp,fs);
    plot(T./60,temp)
    hold on
end
hold off
xlabel('Time (min)')
ylabel('Photobleaching Rate (%)')

% exportgraphics(gcf,'timetraces_asap3_445nm_zoom.pdf', 'ContentType','vector');

%% Load ePhys data

disp('Loading ePhys data...');

fileInfo=dir(fullfile(path,'*.rhd'));

eData=[];eTTL=[];accel=[];

for iFile=1:numel(fileInfo)
    
    read_Intan_RHD2000_file(fullfile(fileInfo(iFile).folder,fileInfo(iFile).name));
    
    eData=[eData; amplifier_data'];
    eTTL=[eTTL; board_dig_in_data'];
    
    if exist('aux_input_data')
        % resampling accelerometer data as fs_accel=eFs/4 by design.
        temp = resample(aux_input_data',4,1);
        accel=[accel; temp];
    end
    
end

eFs=frequency_parameters.board_adc_sample_rate;
eData=notchFilter(eData,eFs,60,'harmonics',4);


%% Data modality alignment

[eData_corr, eTTL_aligned] = ...
    alignMultimodalClockDrift( ...
        eTTL(:,1),...
        oTTL(:,1), ...
        eData,eFs,fs);
    
id_eTTL=find(diff(eTTL_aligned)>0);%id_eTTL=id_eTTL(1:end-1); 
id_oTTL=find(diff(oTTL(:,1))>0);
temp=(id_eTTL-id_oTTL)/fs;

T=id_eTTL;

figure;%
subplot(221)
plot(T/fs,temp*1000,'+b')
axis tight
ylabel('Time (ms)')

subplot(222)
histogram(temp*1000,10)
axis tight
xlim([-3 3])
xlabel('Time (ms)')

subplot(2,2,[3 4])
plot(getTime(eTTL_aligned,fs),[oTTL(:,1) eTTL_aligned]+[0 -1]);
axis tight
% xlim([-3 3])
xlabel('Time (s)')

T=getTime(eData_corr,fs);
figure; plot(T(1:2:end)./60,eData_corr(1:2:end,:))

%% Continue analysis over uncorrupted data segment


if mouseN==1
    % mouse 1: recording glitch: had to stop after ~10 min
    rng=100*fs:size(signal,1);
elseif mouseN==2
    % mouse 2: recording glitch at 600 sec
    rng=600*fs:size(signal,1);
else
    rng=30*fs:size(signal,1); % illumination turned on >20s
end

tpo=signal(rng,:);
lfp=eData_corr(rng,:);

T=getTime(lfp,fs);
figure; plot(T(1:2:end)./60,zscore([tpo(1:2:end,:) lfp(1:2:end,:)])+5.*[0:4])

%% oPhys unmixing
tpo=notchFilter(tpo,fs,297.6,'harmonics',2);
tpo=bpFilter1D(tpo,[inf 300],fs,'order',4);
tpo=runPhotoBleachingRemoval(tpo,'lpCutOff',0.1,'samplingRate',fs,'filterOrder',1)';

plotPSD(tpo,'samplingRate',fs,'BW',[0.5 50],'window',5,'scaleAxis','linear');

figure;plot(getTime(tpo,fs)/fs,bpFilter1D(tpo,[0.5 5],fs))

sig=tpo(:,1);
ref=tpo(:,3);

[umx]=umxCONV(sig,ref,fs,'epoch',5);

M=-100*[umx ref];

plotPSD(lfp,'samplingRate',fs,'BW',[0.5 50],'window',10);
plotPSD(M,'samplingRate',fs,'BW',[0.5 50],'window',10);

fig=figure('Name','LFP Coherence','DefaultAxesFontSize',14,'color','w');
plot(getTime(M,fs),zscore(bpFilter1D([M(:,1) lfp(:,1)],[0.5 25],fs))+[0 5])
xlabel('Time (s)')
ylabel('zscore')

% exportgraphics(gcf,'timetraces.pdf', 'ContentType','vector');

%% Compute LFP coherence plots

idxLFP=2; % pick between channel 1 and 2

n=numel(M(:,1));
idShu=randperm(n);shu=M(idShu,1); %generate the shuffle indices

Csig=[];Csh=[];Cref=[];rng=[];err='sem';

nEpoch=30; win=5; xl=[0.1 100];% set the coherence bandwidth
id=round(linspace(0,n,nEpoch)); disp(['mean epoch duration: ' num2str(mean(diff(id/fs))) ' s'])

% compute frequency vector
rng=id(1)+1:id(1+1);
[~,~,F]=coherencePair(lfp(rng,idxLFP),M(rng,1),fs,'window',win,'plot',false);

parfor i=1:nEpoch-1
    rng=id(i)+1:id(i+1);
    [Csig(:,i)]=coherencePair(lfp(rng,idxLFP),M(rng,1),fs,'window',win,'plot',false);
    [Csh(:,i)]=coherencePair(lfp(rng,idxLFP),shu(rng),fs,'window',win,'plot',false);
    [Cref(:,i)]=coherencePair(lfp(rng,idxLFP),M(rng,2),fs,'window',win,'plot',false);
end

stat=[];

fig=figure('DefaultAxesFontSize',14,'color','w');
plotErrorBar1(Cref,'x_axis',F,'error',err,'figHandle',fig,...
    'color_line',geviColor('ref'),'color_area',geviColor('ref'));
hold on
plotErrorBar1(Csh,'x_axis',F,'error',err,'figHandle',fig,...
    'color_line',geviColor('pace'),'color_area',geviColor('pace'));
hold on
plotErrorBar1(Csig,'x_axis',F,'error',err,'figHandle',fig,...
    'color_line',geviColor('ace'),'color_area',geviColor('ace'));
hold on

id100=find(F==200);
for i=1:id100
    stat(i)=signrank(Csig(i,:),Cref(i,:)); % paired test - observations are paired
end

hold on
alpha=0.05;
stat(stat>alpha)=nan;stat(stat<alpha)=1;
plot(F(1:id100),0*stat,'|k')
hold off
ylabel('Coherence')
xlabel('Time (s)')
legend('ref','','shu','','pv')
xlim(xl);

% exportgraphics(gcf,'coherence_2sWin.pdf', 'ContentType','vector');

%% spindle detection

% M=-100.*M;
n=numel(M(:,1));
idShu=randperm(n);shu=M(idShu,1); %generate the shuffle indices

LFP=lfp(:,1);

% get spindle power
[sw,f] =cwt(LFP,fs,'FrequencyLimits',[1 100],'amor');

id_on=find(f>10,1,'first');
id_off=find(f<15,1,'last');

spinPW_Trace=rescale(mean(abs(sw(id_on:id_off,:))),0,1);
%%
[allSpindles,ttlSpindle]=findSpindles(LFP,'minPeakProm',30);

baselinePrePost=0.5;
stimLength=0.01;

output_swrPow=rasterERP(spinPW_Trace',ttlSpindle,fs,'baselinePrePost',baselinePrePost,'stimLength',stimLength,'figure',false);

%% Parse spindle-related events

plotBW=[0.5 100];
output_LFP=rasterERP(bpFilter1D(LFP,plotBW,fs),ttlSpindle,fs,'baselinePrePost',baselinePrePost,'stimLength',stimLength,'figure',false);

output_umx=rasterERP(bpFilter1D(M(:,1),plotBW,fs),ttlSpindle,fs,'baselinePrePost',baselinePrePost,'stimLength',stimLength,'figure',false);%caxis([-0.006 0.006])
output_ref=rasterERP(bpFilter1D(M(:,2),plotBW,fs),ttlSpindle,fs,'baselinePrePost',baselinePrePost,'stimLength',stimLength,'figure',false);%caxis([-0.005 0.005])
output_shu=rasterERP(bpFilter1D(shu,plotBW,fs),ttlSpindle,fs,'baselinePrePost',baselinePrePost,'stimLength',stimLength,'figure',false);%caxis([-0.005 0.005])

%% spindle trigger average + statistics

stat=[];

time=getTime(output_LFP.arrayRaw(:,1),fs)-baselinePrePost;

fig=figure('Name','SWR trigger average','DefaultAxesFontSize',16,'color','w');
subplot(211)
plotErrorBar1(output_LFP.arrayRaw,'x_axis',time,'figHandle',fig,'color_area',0.*[1 1 1],'color_line',0.*[1 1 1],'error','sem');
hold on
plot([0 0],[-100 100],'--k')
hold off
axis tight
ylabel('V (uV)')
legend('','spindle')
xlim([-baselinePrePost baselinePrePost])

subplot(212)
plotErrorBar1(5*output_swrPow.arrayRaw+0.1,'x_axis',time,'figHandle',fig,'color_area',0.*[1 1 1],'color_line',0.*[1 1 1]);
hold on
plotErrorBar1(output_ref.arrayRaw,'x_axis',time,'figHandle',fig,'color_area',geviColor('ref'),'color_line',geviColor('ref'),'error','sem');
hold on
plotErrorBar1(output_shu.arrayRaw,'x_axis',time,'figHandle',fig,'color_area',geviColor('pace'),'color_line',geviColor('pace'),'error','sem');
hold on
plotErrorBar1(output_umx.arrayRaw,'x_axis',time,'figHandle',fig,'color_area',geviColor('ace'),'color_line',geviColor('ace'),'error','sem');
hold on
plot([0 0],[-1 1],'--k')
hold on

% apply stats PV vs shuffle
for i=1:numel(time)
    stat(i)=signrank(output_umx.arrayRaw(i,:),output_shu.arrayRaw(i,:),'tail','both'); % paired test - observations are paired
end

alpha=0.01;
stat(stat>alpha)=nan;stat(stat<alpha)=1;
plot(time,1*stat,'*k')

hold off
ylabel('-dF/F (%)')
xlabel('Time (s)')
legend('spindles power','','ref','','shuffle','','pv')
axis tight
xlim([-baselinePrePost baselinePrePost])

%% to plot individual event, one-by-one

fig=figure('Name','LFP Coherence','DefaultAxesFontSize',14,'color','w');
T=getTime(output_umx.arrayRaw(:,1),fs);
xlabel('Time (min)')
ylabel('zscore')

%%
i=11; %increment here

plot(time,zscore(bpFilter1D(output_LFP.arrayRaw(:,i),[1 30],fs)))
hold on
plot(time,zscore(bpFilter1D(output_umx.arrayRaw(:,i),[1 5],fs)))
hold off
axis tight









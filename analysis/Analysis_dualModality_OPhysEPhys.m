%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% ANALYSIS of mPFC PV-INs uSMAART/LFP data %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all

% pick mouse number 1 or 2
mouseN=1;

path=['yourPath to uSMAART .mat & LFP .rhd data']; % !!!change this line!!!

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

T=getTime(allData,fs);
figure; plot(T(1:100:end)./60,allData(1:100:end,1:5))

%%

if mouseN==1
% mouse 1: recording glitch: had to stop after ~10 min
rng=840*fs:numel(T);
else
% mouse 2: recording glitch at 600 sec
rng=600*fs:numel(T);
end

oTTL=allData(rng,6:7);
signal=allData(rng,3:5);

T=getTime(signal,fs);
figure; plot(T(1:2:end)./60,signal(1:2:end,:))

%% plot photobleaching time course
figure('DefaultAxesFontSize',14,'color','w');
for i=1:3
temp=signal(:,i);temp=(temp)./median(temp(1:2*fs));
T=getTime(temp,fs);
plot(T./60,temp)
hold on
end
hold off
xlabel('Time (min)')
ylabel('Photobleaching Rate (%)')% exportgraphics(gcf,'timetraces_asap3_445nm_zoom.pdf', 'ContentType','vector');

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
eData=sh_NotchFilter(eData,eFs,60,'harmonics',4);

%% Data modality alignment

[range_ttlE,range_ttlO]=alignSyncTTL(eTTL(:,1),oTTL(:,1),fs,'plotFigure',true);

eTTL_new=eTTL(range_ttlE,:);
eData_new=eData(range_ttlE,:);
if exist('aux_input_data')
    accel=accel(range_ttlE,:);
end

oTTL_new=oTTL(range_ttlO,:);
signal_new=signal(range_ttlO,:);

T=find(diff(eTTL_new(:,1))>0);
temp=(find(diff(eTTL_new(:,1))>0)-find(diff(oTTL_new(:,1))>0))/fs;

figure;
plot(T/fs,temp,'+b')
axis tight

lfp=eData_new;
plotPSD(lfp,'samplingRate',fs,'BW',[0.5 50],'window',5,'scaleAxis','linear');
figure;plot(getTime(lfp,fs),bpFilter1D(lfp,[0.5 5],fs))

%% oPhys unmixing
signal_new=sh_NotchFilter(signal_new,fs,297.6,'harmonics',2);
signal_new=bpFilter1D(signal_new,[inf 300],fs,'order',4);
signal_new=runPhotoBleachingRemoval(signal_new,'lpCutOff',0.1,'samplingRate',fs,'filterOrder',1)';

plotPSD(signal_new,'samplingRate',fs,'BW',[0.5 50],'window',5,'scaleAxis','linear');

figure;plot(getTime(signal_new,fs)/fs,bpFilter1D(signal_new,[0.5 5],fs))

sig=signal_new(:,1);
ref=signal_new(:,3);

[umx]=umxCONV(sig,ref,fs,'epoch',5);

%% 

if mouseN==1
    % mouse 1
    rng=1:numel(lfp(:,1));
else
    % mouse 2
    rng=100*fs:1900*fs;
end

lfp=lfp(rng,:);
M=[umx ref];M=M(rng,:);

plotPSD(lfp,'samplingRate',fs,'BW',[0.5 50],'window',10);
plotPSD(M,'samplingRate',fs,'BW',[0.5 50],'window',10);

fig=figure('Name','LFP Coherence','DefaultAxesFontSize',14,'color','w');
plot(getTime(M,fs),zscore(bpFilter1D([M(:,1) lfp(:,1)],[0.5 25],fs)))
xlabel('Time (min)')
ylabel('zscore')

% exportgraphics(gcf,'timetraces.pdf', 'ContentType','vector');

%% LFP - PV - Ref coherence

fig=figure('Name','LFP Coherence','DefaultAxesFontSize',14,'color','w');
coherencePair(lfp(rng,1),umx(rng),fs,'window',5,'figHandle',fig);
hold on
coherencePair(lfp(rng,1),ref(rng),fs,'window',5,'figHandle',fig);
hold off
xlim([0.1 50])
title('Coherence LFP #1')

fig=figure('Name','LFP Coherence','DefaultAxesFontSize',14,'color','w');
coherencePair(lfp(rng,2),umx(rng),fs,'window',5,'figHandle',fig);
hold on
coherencePair(lfp(rng,2),ref(rng),fs,'window',5,'figHandle',fig);
hold off
xlim([0.1 50])
title('Coherence LFP #2')


%% Compute LFP coherence plots

idxLFP=1; % pick between channel 1 and 2

n=numel(M(:,1));
idShu=randperm(n);shu=M(idShu,1); %generate the shuffle indices

Csig=[];Csh=[];Cref=[];rng=[];err='sem';

nEpoch=30; win=2; xl=[0.1 100];% set the coherence bandwidth
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









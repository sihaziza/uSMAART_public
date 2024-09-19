
%% Prepare all path and gather ephys and 1 oPhys data

mouse='m915';
dateMeas='20210221';
rawEPhysPath='F:\GEVI_Spike\ePhys\Spontaneous';
rawOPhysPath='F:\GEVI_Spike\Preprocessed\Spontaneous';
savePath='F:\GEVI_Spike\Analysis\Spontaneous';

if ~exist(fullfile(savePath,mouse,dateMeas),'dir')
    mkdir(fullfile(savePath,mouse,dateMeas));
end

% Load ePhys data
ePhysPath=fullfile(rawEPhysPath,mouse,dateMeas);%'F:\GEVI_Spike\ePhys\Spontaneous\m915\20210221';

ePhysPath=dir(fullfile(ePhysPath,'*.rhd'));
fprintf('%2.0f intan files detected \n',length(ePhysPath));
disp('Recursive readout of INTAN data')

nEphysFile=1;
read_Intan_RHD2000_file(fullfile(ePhysPath(nEphysFile).folder,ePhysPath(nEphysFile).name))
ePhys.data=amplifier_data';
ePhys.ttl=board_dig_in_data';
ePhys.fps=frequency_parameters.board_adc_sample_rate;

% Find all oPhys measurements folder
oPhysPath=fullfile(rawOPhysPath,mouse,dateMeas);%'F:\GEVI_Spike\Preprocessed\Spontaneous\m915\20210221';
% oPhysPath='F:\GEVI_Spike\Preprocessed\Spontaneous\m913\20210215';
folder=dir(fullfile(oPhysPath,'meas*'));

% look at the first oPhys measurements, assuming all meas are identical

% folder=strcat('meas0',num2str(iMeas-1));
filePath=fullfile(oPhysPath,folder(end).name,'metadata.mat');
% Load Behavior data
if isfile(filePath)
    load(filePath);
    oPhys.data=metadata.Locomotion; %> change with single Units~
    oPhys.ttl=metadata.TTL;
    oPhys.loco=metadata.Locomotion;
    oPhys.fps=metadata.fps;
    [~,oPhys.speed,oPhys.locoRestTTL,~]=getMouseSpeed(oPhys.loco,oPhys.fps);
else
    warning('No behavioral data - is it normal?')
end

% trace=zscore(ica_sig(end,:))';
% theta=sh_bpFilter(trace,[0.1 100],oPhys.fps);
% z=hilbert(theta);
% phase=angle(z);
% temp=-[diff(phase); 0];
% [pks,locs] =findpeaks(temp,'MinPeakHeight',pi);
% for i=1:numel(pks)
%     try
%    sortedPhase(:,i)= trace(locs(i)-round(fs/3):locs(i)+round(fs/3));
%     end
% end
% plot(sortedPhase)
% plot([phase temp])
% ai=find(diff(phase)>pi/4);
% plot(
% af=find(phase==max(phase));
% t=getTime(oPhys.locoRestTTL,oPhys.fps);
% fs= oPhys.fps;
% instF = fs/(2*pi)*diff(unwrap(angle(z)));
% instF2=instfreq(theta,t,'Method','hilbert');figure();plot(instF2)
% plot(getTime(instF,fs),instF)
% plot(t,[oPhys.locoRestTTL oPhys.speed trace-4 angle(z)-2])
%% Filter ePhys data assuming fpsO<fpsE

ePhysDuration=(length(ePhys.ttl)-1)/ePhys.fps;
timeOld=0:1/ePhys.fps:ePhysDuration; % real fpsE sampling
timeNew=0:1/oPhys.fps:ePhysDuration; % sample ePhys at fpsO

% Filtering at newFPS/2 to fulfill Shannon theorem
temp=sh_bpFilter(ePhys.data,[inf oPhys.fps/2],ePhys.fps);
ePhysDS.data=interp1(timeOld,temp,timeNew);
ePhysDS.ttl=interp1(timeOld,ePhys.ttl,timeNew);
% Correct for 0/1 transition after interpolation
temp=ePhysDS.ttl;
temp=double(temp>0.5);
ePhysDS.ttl=double(~(temp<0.5))';

ePhysDS.fps=oPhys.fps;
fps=oPhys.fps;
h=figure();
plot(timeOld,ePhys.ttl,timeNew,ePhysDS.ttl,'o','markersize',8)
savePDF(h,'ePhys TTL interp',fullfile(savePath,mouse,dateMeas))

% Now time to align ePhys/oPhys
% Assigne ePhys data into chunck corresponding to oPhys recordings
[parsedData]=parseEPhys(ePhysDS,'savePath',fullfile(savePath,mouse,dateMeas));

% Check that the number of parsed data equal number of oPhys measurements
if numel(parsedData)~=numel(folder)
    error('Houston, we have a problem... could not find equal number of measurements...')
else
    disp('Everything looks good to me!')
end

%% If everything looks right, batch process all files.
for iMeas=1:numel(folder)
    folderMeas=strcat('meas0',num2str(iMeas-1));
    filePath=fullfile(oPhysPath,folderMeas,'metadata.mat');
    tempSave=fullfile(savePath,mouse,dateMeas,folderMeas);
    if ~exist(tempSave,'dir')
        mkdir(tempSave);
    end
    
    % Load Behavior data
    if isfile(filePath)
        load(filePath);
        oPhys.data=metadata.Locomotion; %> change with single Units~
        oPhys.ttl=metadata.TTL;
        oPhys.loco=metadata.Locomotion;
        oPhys.fps=metadata.fps;
        [~,oPhys.speed,oPhys.locoRestTTL]=getMouseSpeed(oPhys.loco,oPhys.fps,'savePath',tempSave);
        
    else
        warning('No behavioral data - is it normal?')
    end
    
    % Align syncTTL from oPhys and ePhys
    [range_ttlE,range_ttlO]=alignSyncTTL(parsedData(iMeas).ttl,oPhys.ttl,parsedData(iMeas).fps,'savePath',tempSave);
    
    % Final assignment of the correct data points [output]=alignRecording(ePhysDS,oPhys);
    output.dataE=parsedData(iMeas).data(range_ttlE,:);
    output.ttlE=parsedData(iMeas).ttl(range_ttlE);
    output.dataO=oPhys.data(range_ttlO);
    output.ttlO=oPhys.ttl(range_ttlO);
    output.speed=oPhys.speed(range_ttlO);
    output.brainStates=oPhys.locoRestTTL(range_ttlO);
    output.time=getTime(output.ttlO,fps)';
    output.fps=fps;
    
    outputName=['output_alignedEPhysOPhys_' datestr(datetime('now'),'yyyymmddTHHMMSS') '.m'];
    save(fullfile(tempSave,outputName),'output');
    close all;
end

%%
% velocity=rescale(output.speed
imagesc(M')
caxis([0 5])
colormap(redblue)
axis off;

%%
lfp=output.dataE(:,2);

[wt,f] = cwt(lfp,output.fps);

%%

BW=[5 10];
% descending order here...
[highBWidx,~,~]=find(f>=BW(1),1,'last');
[lowBWidx,~,~]=find(f<=BW(2),1,'first');
% output is complex > take the absolute value. Quid phase?
power=abs(wt);
powerDelta=mean(power(lowBWidx:highBWidx,:),1)';

BW=[4 8];
% descending order here...
[highBWidx,~,~]=find(f>=BW(1),1,'last');
[lowBWidx,~,~]=find(f<=BW(2),1,'first');
% output is complex > take the absolute value. Quid phase?
power=abs(wt);
powerTheta=mean(power(lowBWidx:highBWidx,:),1)';

plot(getTime(powerTheta,output.fps),[output.speed powerDelta powerTheta])
%%

Mtheta=[powerTheta powerTheta];
Mloco=[output.speed output.speed];

figure()
subplot(211)
imagesc(Mloco')
caxis([0 5])
colormap(parula)
axis off;

subplot(212)
imagesc(Mtheta')
% caxis([50 80])
colormap(parula)
axis off;

%%

%%
M=[signal reference umx lfp1 lfp2];

Mfilt=sh_bpFilter(M,[0.5 50],fs);

plotPSD(M,'FrameRate',outputParse.fs,'FreqBand',[0.1 150],'Window',5);
legend('Green','Red','umx','lfp1','lfp2')

plot(getTime(ttl,fs),[Mfilt ttl]);
[powerTrace]=plotPowerEvolution(umx,fs);

%%
% Mfilt=sh_bpFilter(M,[0.5 50],fs);
[outputParse]=getStimEpoch(umx,ttl,fs,'baselinePrePost',1,'getShuffle',true,'preStimNorm',true);

% Make sliding window bandpower!

epoch=outputParse.stimBand;
wave(:,:,1)=outputParse.arrayRaw;
wave(:,:,2)=outputParse.arrayShuffle;
%%
nn=size(wave,3);
% nameChannel={'Green','Red','Umx','lfp1','lfp2'};
nameChannel={'ERP','Shuffled'};
% find the min/max from ERP, not shuffle
% temp=wave(:,:,1);
% cRange=[min(temp(:)) max(temp(:))];
cRange=[-5 5];

figure(100)
for i=1:nn
    subplot(1,nn,i)
    sh_singleERPcolor(wave(:,:,i),epoch,cRange,fs)
    title(nameChannel{i})
    xlim(epoch)
end

temp=mean(wave(:,:,1),2);
cRange=[min(temp(:)) max(temp(:))];

figure(110)
for i=1:nn
    temp=mean(wave(:,:,i),2);
    
    tm=(epoch(1):1/fs:epoch(2)-1/fs)';
    
    subplot(1,nn,i)
    plot(tm,temp,'k')
    ylim(cRange)
    hold on
    plot([0 0],[min(temp(:)) max(temp(:))],'--k','LineWidth',1)
    hold off
    xlim(epoch)
    
end

stat={'std','c99','c95','sem'};
h=figure('DefaultAxesFontSize',18,'color','w');

for i=1:numel(stat)
    options.handle     = h;
    options.alpha      = 0.5;
    options.line_width = 2;
    options.error      = stat{i};
    
    subplot(1,numel(stat),i)
    options.color_area = [128 193 219]./255;    % Blue theme
    options.color_line = [ 52 148 186]./255;
    plot_areaerrorbar(wave(:,:,1),options)
    hold on
    options.color_area = [243 169 114]./255;    % Orange theme
    options.color_line = [236 112  22]./255;
    plot_areaerrorbar(wave(:,:,2),options)
    hold off
    title(options.error)
end

% get spectrogram then theta, delta power same as speed
%%
% plot speed vs ePhys in the delta and theta range to assess coarse
% correlation
delta=zscore(sh_bpFilter(output.dataE,[1 4],output.fps));
theta=zscore(sh_bpFilter(output.dataE,[6 10],output.fps));

figure()
plot(output.time,[delta theta-5 zscore(output.dataE)-10 output.speed+5 output.brainStates+5])

tic;plotPSD(ePhys.data,'FreqBand',[0.5 100],'FrameRate',ePhys.fps);toc;
tic;plotPSD(ePhysDS.data,'FreqBand',[0.5 100],'FrameRate',ePhysDS.fps);toc;

%%
timeE=getTime(ttlE,fpsE);
timeO=getTime(ttlO,fpsO);

plot(timeE, [ttlE 10*movmean(ttlE,30*fpsE)], timeO,ttlO - linspace(0,11,6))

[~,locs] =findpeaks(movmean(ttlE,30*fpsE),'MinPeakDistance',30*fpsE);
delta=diff(locs)./2;
delta=[0 delta' length(ttlE)-locs(end)];

for iLocs=1:length(locs)
    meas(iLocs).ttlE= ttlE(round(locs(iLocs)-delta(iLocs)):round(locs(iLocs)+delta(iLocs+1)));
end

for iMeas=1:6
    ttlE_temp=meas(iMeas).ttlE;
    ttlO_temp=ttlO(:,iMeas);
    
    x=getTime(ttlO_temp,fpsO);
    v=ttlO_temp;
    xq=linspace(0,x(end),x(end)*fpsE);
    vq = interp1(x,v,xq,'nearest')';
    plot(xq, vq,ttlE_temp)
    
    % Find the delay between the last ttl and the end of the recording
    [offsetO,~]=find(diff(ttlO_temp)==1,1,'last');
    dT=(length(ttlO_temp)-offsetO+1)/fpsO;
    
    % detect the equivalent end on the ePhys data
    [offsetE,~]=find(diff(ttlE_temp)==1,1,'last');
    %     truncE=ttlE_temp(
    timeE=getTime(ttlE_temp,fpsE);
    timeO=getTime(ttlO_temp,fpsO);
    
    [onsetO,~]=find(diff(ttlO_temp)==1,1,'first');
    [offsetO,~]=find(diff(ttlO_temp)==-1,1,'last');
    timeO([onsetO offsetO end])
    
    [onsetE,~]=find(diff(ttlE_temp)==1,1,'first');
    [offsetE,~]=find(diff(ttlE_temp)==-1,1,'last');
    timeE([onsetE offsetE end])
    
    plot(timeE, ttlE_temp, timeO,ttlO_temp)
    [allEventsO,~]=find(diff(ttlO_temp)==1);
    [allEventsE,~]=find(diff(ttlE_temp)==1);
    
    if numel(allEventsO)~=numel(allEventsE)
        warning('n')
    end
end
%%

% Check for incorrect parsing of data
clear ePhys oPhys dio metadata

ePhys=amplifier_data(:,idE_first:idE_first+delta-1)';
oPhys=datatemp(idO_first:idO_first+delta-1,1:4);
dio=datatemp(idO_first:idO_first+delta-1,dioChannelO);

figure()
subplot(311)
plot(time,zscore(ePhys))
title('LFP recordings')

subplot(312)
plot(time,zscore(oPhys))
title('TEMPO recordings')

subplot(313)
plot(time,dio,'k')
title('Camera Frames')

ylim([0 1.1])
xlabel('Time (s)')

% Save Metadata
metadata.fps=min(fpsO,fpsE);
metadata.dioInput='true';
metadata.dioInputChannel={'AirPuff'};
metadata.ePhysChannel={'LFP1' 'LFP2' 'LFP3'};
metadata.oPhysChannel={'DIO-Ace1.0' 'cag-cyOFP' 'Camk2-Varnam1.0'};
metadata.mouseID='m4 PV-Cre Het male';
metadata.mouseDOB='20191107';
metadata.originalFilePath=directory;
metadata.originalePhysFile=ePhysPath.name;
metadata.originaloPhysFile=oPhysPath(2).name;
metadata

% Save parsed data
pathName=fullfile(metadata.originalFilePath,'preprocessed');
if ~exist(pathName,'dir')
    mkdir(pathName)
end
save(fullfile(pathName,'dio.mat'),'dio');
save(fullfile(pathName,'metadata.mat'),'metadata');
save(fullfile(pathName,'ePhys.mat'),'ePhys');
save(fullfile(pathName,'oPhys.mat'),'oPhys');
%%



Here is the code I used for the first movie:





clear;

%M = hdf5read('yale_NDNF-VIP_moco.h5','/movie');

M = hdf5read('m83_d201124_s06dualColorSlidePulsingLEDs-fps781-cR_moco.h5','/mov');







config=[];

config=get_defaults(config);

config.spatial_highpass_spatial=inf;

config.avg_cell_radius=10;

M=preprocess_movie(M,config);

M=abs(M);

X=ones(7,7,1)/49; % moving average of 7x7 square, for the other movie 3x3 was good.

M=convn(M,X,'same');

%M=-M;



%%
directory='F:\GEVI_Spike\Preprocessed\Spontaneous\m915\20210221';
for i=0:4
    meas=['meas0' num2str(i)];
    path=dir(fullfile(directory,meas,'*_dtr.h5'));
    pathG=fullfile(path(1).folder,path(1).name);
    tic;runEXTRACT(pathG,'polarityGEVI','pos','checkCell',false);toc;
    tic;runEXTRACT(pathG,'polarityGEVI','neg','checkCell',false);toc;
    pathR=fullfile(path(2).folder,path(2).name);
    tic;runEXTRACT(pathG,'polarityGEVI','neg','checkCell',false);toc;
end

directory='F:\GEVI_Spike\Preprocessed\Spontaneous\m915\20210221';
for i=5
    meas=['meas0' num2str(i)];
    path=dir(fullfile(directory,meas,'*_dtr.h5'));
    pathG=fullfile(path(1).folder,path(1).name);
    tic;[M]=runCheckCell(pathG,'polarityGEVI','pos','checkCell',false);toc;
%     tic;runEXTRACT(pathG,'polarityGEVI','neg','checkCell',false);toc;
%     pathR=fullfile(path(2).folder,path(2).name);
%     tic;runEXTRACT(pathG,'polarityGEVI','neg','checkCell',false);toc;
end
%%
a='F:\GEVI_Spike\Preprocessed\Visual\m81\20210402\meas00\DemixingEXTRACT\m81_d210402_s00driftgrating315-fps750-cG_bp_moco_dtr\extractOutput_neg_20210402T160704.m';

% a='F:\GEVI_Spike\Preprocessed\Spontaneous\m915\20210221\meas05\DemixingEXTRACT\m915_d210221_s05laser100pct--fps601-cG_bp_moco_dtr';
% b=fullfile(a,'extractOutput_pos_20210401T130938.m');
b_save=strrep(a,'.m','_cleaned.mat');
output=importdata(a);
% cell_check(output,M)
%%
% m=full(spatial_weights(:,:,1));
nUnits=10;
% decade=1;
figure(1)
for i=1:nUnits
subplot(nUnits,nUnits,[nUnits*(i-1)+1 nUnits*(i-1)+2])
imagesc(full(output.spatial_weights(:,:,(decade-1)*nUnits+i)))
ylabel(num2str((decade-1)*nUnits+i))
subplot(nUnits,nUnits,[nUnits*(i-1)+3 nUnits*i])
plot(output.temporal_weights(:,(decade-1)*nUnits+i))
end
decade=decade+1;
%%
% list=[1 2 3 4 8 9 10 11 13 14 20 31 39 41 91 137];
% list_down=[52 55 58 66 97 120];
list=[1 2 5 7 8 13 17 20 31 32];
nUnits=numel(list);
sFilt=zeros(size(output.spatial_weights,1),size(output.spatial_weights,2),nUnits);
tFilt=zeros(size(output.temporal_weights,1),nUnits);
figure(1)
for i=1:nUnits
% subplot(nUnits,nUnits,[nUnits*(i-1)+1 nUnits*(i-1)+2])
sFilt(:,:,i)=rescale(full(output.spatial_weights(:,:,list(i))),0,10);
tFilt(:,i)=output.temporal_weights(:,list(i));

% imagesc(full(spatial_weights(:,:,list(i))))
% ylabel(num2str(list(i)))
% subplot(nUnits,nUnits,[nUnits*(i-1)+3 nUnits*i])
%     plot(temporal_weights(:,list(i)))
end
figure(2)
imagesc(sum(sFilt,3))

%%
% save fps as well...
neg=load('F:\GEVI_Spike\Preprocessed\Spontaneous\m915\20210221\meas05\DemixingEXTRACT\m915_d210221_s05laser100pct--fps601-cG_bp_moco_dtr\extractOutput_neg_20210401T124445_cleaned.mat');
pos=load('F:\GEVI_Spike\Preprocessed\Spontaneous\m915\20210221\meas05\DemixingEXTRACT\m915_d210221_s05laser100pct--fps601-cG_bp_moco_dtr\extractOutput_pos_20210401T130938_cleaned.mat');

time=getTime(pos.output.temporal_filter(:,1),601);
 %%
 tot=size(pos.output.temporal_filter,2)+size(neg.output.temporal_filter,2);
 steps=linspace(0,0.04*tot,tot);
 M=[];
figure(1)
for i=1:size(pos.output.temporal_filter,2)
M(:,i)=rescale(pos.output.temporal_filter(:,i),0,0.05);
end
plot(time,M+steps(1:size(M,2)),'k')
hold on
 M=[];

for i=1:size(neg.output.temporal_filter,2)
M(:,i)=rescale(neg.output.temporal_filter(:,i),0,0.05);
end
plot(time,M+steps(end-size(M,2)+1:end),'r')
hold off
disp('hello')
xlim([-1 31])
ylim([-0.1 0.9])
%%
clear output
output.spatial_filter=sFilt;
output.temporal_filter=tFilt;
output.original_listID_pos=list;
output.original_listID_neg=list;

save(b_save,'output');

% figure()
% imagesc(sum(sFilt,3))
% c1=colormap(flipud(cool));c1(1:32,:)=ones(32,3);
% c2=colormap((cool));c2(1:32,:)=ones(32,3);
% colormap(parula)
%%
%load movie here

config=[];

config = get_defaults(config); %calls the defaults

config.dendrite_aware=1; %If you want dendrites from the movie

config.use_gpu=1;



% Essentials, without these EXTRACT will give an error:

config.avg_cell_radius=5; %Pick a reasonable cell radius





%Optionals, but strongly advised to handpick:

%Movie is small enough that EXTRACT will not automatically partition this,

%but still a good idea to keep these in sight!

config.num_partitions_x=1;

config.num_partitions_y=1;

config.preprocess=1;

config.spatial_highpass_spatial=inf;

config.cellfind_filter_type='none'; % Feel free to use your own filters,

%but I do suggest trying 'butter' at times.



config.max_iter=6;



config.cellfind_max_steps=50;



% Voltage specific configs

config.trace_output_option='nonneg';

config.kappa_std_ratio=100;

config.cellfind_kappa_std_ratio=100;

config.cellfind_min_snr=0;

config.thresholds.T_min_snr=0;



% Perform EXTRACTion:

output=extractor(M,config);








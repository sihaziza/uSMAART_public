
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% analyze time jitter of visual stimulation %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all

% set the data path
mainPath='Y:\Hansol\Opioid';
folder='WT_delta_Rcamp_BLA3215_autosave_000';
path=fullfile(mainPath,folder);
file=dir(fullfile(path,'*.mat'));

data=[];

% load uSMAART data
nbytes = fprintf('processing %d%% \n', 0);
nFile=numel(file);
for iFile=1:nFile
    [temp_pre,meta]=loaduSMAART2mat(fullfile(file(iFile).folder,file(iFile).name),'figure',false,'verbose',false);
    data=[data; temp_pre];
    
    fprintf(repmat('\b',1,nbytes))
    nbytes = fprintf('processing %d%% \n', round(iFile/nFile*100));
end

[u,v]=find(isnan(data));
data=data(setdiff(1:length(data),unique(u)),:);

fs=meta.fs;

start=30000;stop=size(data,1);
ttl=data(start:stop,6);

figure;
plot(getTime(data,fs), data(:,3:7))
legend('green','cyofp','red','airpuf','tms')

%% unmixing

signals=data(start:stop,3:5);

% data conditioning
signals=sh_NotchFilter(signals,fs,297.6,'harmonics',2);
signals=bpFilter1D(double(signals),[inf 200],fs,'order',2);
signals=runPhotoBleachingRemoval(signals,'lpCutOff',0.1,'samplingRate',fs,'filterOrder',2)';

sig=signals(:,1); 
ref=signals(:,2); % ref channel is either #2 or #3 for double or single sensor, resp.

% unmxing step
epoch=1.5;
opi=umxCONV(sig, ref, fs,'epoch',epoch);

% plot PSD pre vs post unmixing
plotPSD(signals,'samplingRate',fs,'BW',[1 100],'window',5);
plotPSD([sig ref opi],'samplingRate',fs,'BW',[1 100],'window',5);

%% plot event-triggered rasters and average

baselinePrePost=1;
stimLength=2;

[output_sig]=rasterERP(sig*100,ttl,fs,'baselinePrePost',baselinePrePost,'stimLength',stimLength,'figure',true);
[output_ref]=rasterERP(ref*100,ttl,fs,'baselinePrePost',baselinePrePost,'stimLength',stimLength,'figure',true);
[output_opi]=rasterERP(opi*100,ttl,fs,'baselinePrePost',baselinePrePost,'stimLength',stimLength,'figure',true);

fig=figure('Name','Trigger average','DefaultAxesFontSize',16,'color','w');
time=getTime(output_opi.arrayRaw(:,1),fs)-baselinePrePost;
plotErrorBar1(output_sig.arrayRaw,'x_axis',time,'figHandle',fig,'color_area',geviColor('ace'),'color_line',geviColor('ace'),'error','c95');
hold on
plotErrorBar1(output_ref.arrayRaw,'x_axis',time,'figHandle',fig,'color_area',geviColor('vnm'),'color_line',geviColor('vnm'),'error','c95');
hold on
plotErrorBar1(output_opi.arrayRaw,'x_axis',time,'figHandle',fig,'color_area',geviColor('pace'),'color_line',geviColor('pace'),'error','c95');
hold on
plot([0 0],[-0.5 1.5],'--k')
hold off
ylabel('dF/F (%)')
xlabel('Time (s)')
% legend('red','','ref','','rcamp')
legend('green','','ref','','opioid')

%%


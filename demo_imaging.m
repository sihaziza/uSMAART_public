clear all; close all;

mainPath='Y:\GEVI_Wave\Preprocessed\Visual';

mousePath='m883038\20250327\meas00';

fileInfo=dir(fullfile(mainPath,mousePath,'*.h5'));

gName=fullfile(fileInfo(2).folder,fileInfo(2).name);fileInfo(2).name
rName=fullfile(fileInfo(5).folder,fileInfo(5).name);fileInfo(5).name

fs=196;

savePath=['H:\Analysis_endothelial_20240131\' mousePath];

if ~isfolder(savePath)
   mkdir(savePath) 
end

% 
[~, tG,boundbox]=getUnitsROI(gName,fs);
[~, tR]=getUnitsROI(rName,fs,'boundbox',boundbox);

% load behavior data
fileInfo=dir(fullfile(mainPath,mousePath,'**/*_framestamps 0.txt'));
beha=importdata(fullfile(fileInfo.folder,fileInfo.name));
ttl=beha.data(1:length(tG),4);

%% DETRENDING STEP

tG_dtr=runPhotoBleachingRemoval(tG,'lpCutOff',0.1,'samplingRate',fs,'filterOrder',2')';
tR_dtr=runPhotoBleachingRemoval(tR,'lpCutOff',0.1,'samplingRate',fs,'filterOrder',2)';

Mraw=[tG tR];
figure;plot(getTime(tG,fs),Mraw)

Mraw_dtr=[tG_dtr tR_dtr];
figure;plot(getTime(tG,fs),Mraw_dtr)
plotPSD(Mraw_dtr,'samplingRate',fs,'BW',[1 30],'window',3);

%% UNMIXING STEP

umx=[];
coef=5:5:100;
for i=1:numel(coef)
[umx(:,i)]=umxCONV(tG_dtr,tR_dtr,fs,'epoch',0.1*coef(i)+0.5);
end

% ASAP3 with negative polarity
umx=-(umx-mean(umx)); 

plotPSD([tG_dtr umx],'samplingRate',fs,'BW',[0.01 10],'window',10,'scaleAxis','log2');
legend()

%% After selecing the best unmixing coefficient (defined as the coefficient that minimize aberrant frequency power)

iBest=19;
[umx]=umxCONV(tG_dtr,tR_dtr,fs,'epoch',0.1*coef(iBest-1)+0.5);

Mumx=100.*[tG_dtr tR_dtr umx];
plotPSD(Mumx,'samplingRate',fs,'BW',[0.01 50],'window',10);

fig=figure('DefaultAxesFontSize',14,'color','w');
for i=1:3
    temp=Mumx(:,i);
    [output{i}]=rasterERP(temp,ttl,fs,'baselinePrePost',2,'stimLength',4,'figHandle',fig,'figure',false);
end

fig=figure('DefaultAxesFontSize',14,'color','w');
t=getTime(output{1}.arrayRaw,fs)+output{1}.stimBand(1);
plotErrorBar1(output{2}.arrayRaw,'x_axis',t,'error','sem','figHandle',fig,...
    'color_line',geviColor('ref'),'color_area',geviColor('ref'));
hold on
plotErrorBar1(output{3}.arrayRaw,'x_axis',t,'error','sem','figHandle',fig,...
    'color_line',geviColor('ace'),'color_area',geviColor('ace'));
% hold on
% plotErrorBar1(output{1}.arrayRaw,'x_axis',t,'error','sem','figHandle',fig,...
%     'color_line',geviColor('pace'),'color_area',geviColor('pace'));
hold on
plot([0 0],[-1 2],'--k','LineWidth',1)
hold off
ylabel('{\Delta}F/F (%)')
xlabel('Time (s)')
axis tight
legend('ref','','endo','','asap3')

%% SAVING DATA for further analysis

% data.tG=tG;
% data.tR=tR;
% data.umx=umx;
% data.ttl=ttl;
% data.fs=fs;

% save(fullfile(savePath,'data.mat'),'data');

%%



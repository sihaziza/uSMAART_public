function [restBoot,runBoot,stateTTL]=getRestRunBootIndex(speed,varargin)

%% DEFAULT OPTIONS
options.speedCutOff=2; % in cm/s
options.transCutOff=0.2; % mean speed value to assign 'resting state'
options.maxNumChanges=15;
options.resample=true;
options.samplingRate=2000;

%% UPDATE OPTIONS
if nargin>1
    options=getOptions(options,varargin);
end

%%

[stateTTL]=findRestRunTransitions(speed,...
    'speedCutOff',options.speedCutOff,...
    'maxNumChanges',options.maxNumChanges,...
    'resample',options.resample);

% identify the bimodal state indexes.

id_restRun=find(diff(stateTTL)==1);
id_runRest=find(diff(stateTTL)==-1);

if (stateTTL(1)==0)&&(stateTTL(end)==1)
    restBoot=[[1 id_runRest+1];id_restRun]';
    runBoot=[id_restRun+1;[id_runRest numel(stateTTL)]]';
    
elseif (stateTTL(1)==1)&&(stateTTL(end)==1)
    restBoot=[id_runRest+1;id_restRun]';
    runBoot=[[1 id_restRun+1];[id_runRest numel(stateTTL)]]';
    
elseif (stateTTL(1)==0)&&(stateTTL(end)==0)
    restBoot=[[1 id_runRest+1];[id_restRun numel(stateTTL)]]';
    runBoot=[id_restRun+1;id_runRest]';
    
elseif (stateTTL(1)==1)&&(stateTTL(end)==0)
    restBoot=[id_runRest+1; [id_restRun numel(stateTTL)]]';
    runBoot=[[1 id_restRun+1];id_runRest]';
end

T=getTime(speed,options.samplingRate);

fig=figure('DefaultAxesFontSize',16,'color','w');
plot(T,speed)
hold on
plot(T,stateTTL-1.5) % SH-20221201-changed speed' to speed
hold off
ylabel('Speed (cm/s)')
axis tight
ylim([-2 15])
% title(strrep(name,'_','-'))

% [savePath,name]=fileparts(fileName);
% saveName=fullfile(savePath,[name '.png']);
% saveas(fig,saveName);
%     saveName=fullfile(savePath,[name '.pdf']);
%     exportgraphics(fig,saveName,'ContentType','vector','BackgroundColor','none');


end
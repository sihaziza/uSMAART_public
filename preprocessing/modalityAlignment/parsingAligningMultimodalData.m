function [output]=parsingAligningMultimodalData(path)
% this function is parsing ePhys, oPhys and behavioral data, interpolate
% the data so as to get a common sampling rate and aligning them to a
% common synchronizing TTL

%% PULLING ALL DATA

% Parsing Electrophysiological data
ePhysDataPath=dir(fullfile(path,'*.rhd'));
if ~isempty(ePhysDataPath)
    ePhysData=[];t_E=[];ttl_E=[];
    for iRHD=1:length(ePhysDataPath)
        [output_ePhys]=read_Intan_RHD2000_file(fullfile(ePhysDataPath(iRHD).folder,ePhysDataPath(iRHD).name));
        ePhysData=[ePhysData; output_ePhys.dataAI];
        t_E=[t_E; output_ePhys.timeStamps];
        ttl_E=[ttl_E; output_ePhys.dataDI];
        fs_E=output_ePhys.samplingRate;
    end
else
    disp('No INTAN data found...')
end

% Parsing Fiber Photometry data
oPhysDataPath=dir(fullfile(path,'*TEMPO.mat'));
if ~isempty(oPhysDataPath)
    load(fullfile(oPhysDataPath.folder,oPhysDataPath.name));
    oPhysData=Data(:,1:end-1);
    t_O=TimeStamps;
    ttl_O=Data(:,end);
    fs_O=1/mean(diff(TimeStamps));
else
    disp('No FiberPhotometry data found...')
    oPhysData=[];
end
% Parsing Dehavioral data
behaDataPath=dir(fullfile(path,'*.txt'));
if ~isempty(behaDataPath)
    beha=importdata(fullfile(behaDataPath.folder,behaDataPath.name));
    
    % interpolate behavioral data which has non-fix sampling rate
    ti=beha.data(:,1)/1000;
    loco=beha.data(:,2);
    ttl=beha.data(:,3);
    % remove non-unique points
    [a,~]=find(~diff(ti));
    ti(a)=[];loco(a)=[];ttl(a)=[];
    
    t_B=(0:1/fs_O:ti(end))';
    loco=interp1(ti,loco,t_B,'nearest');
    ttl_B=interp1(ti,ttl,t_B,'nearest');
    fs_B=1/mean(diff(t_B));
else
    disp('No Behavioral data found...')
    loco=[];
end

%% aligning all modalities together
% if all ePhys, oPhys and Behavior data found...
if ~isempty(ePhysDataPath)&&~isempty(oPhysDataPath)&&~isempty(behaDataPath)
    figure(2)
    subplot(211)
    plot(t_B,ttl_B+1,t_E,ttl_E,t_O,ttl_O-1);
    legend('behavior','ePhys','oPhys')
    title('synchronizing TTL pre-alignment')
    
    % aligning ePhys & oPhys
    [range_ttlE,range_ttlO]=alignSyncTTL(ttl_E,ttl_O,fs_E);
    ttl_E=ttl_E(range_ttlE);ePhysData=ePhysData(range_ttlE,:);
    ttl_O=ttl_O(range_ttlO,:);oPhysData=oPhysData(range_ttlO,:);
    
    % aligning ePhys & Behavior, correct oPhys
    [range_ttlE,range_ttlB]=alignSyncTTL(ttl_E,ttl_B,fs_E);
    ttl_E=ttl_E(range_ttlE);
    ePhysData=ePhysData(range_ttlE,:);
    oPhysData=oPhysData(range_ttlE,:);ttl_O=ttl_O(range_ttlE,:);
    loco=loco(range_ttlB);ttl_B=ttl_B(range_ttlB);
    
    subplot(212)
    plot(getTime(ttl_E,fs_E),[ttl_B+1 ttl_E ttl_O-1]);
    legend('behavior','ePhys','oPhys')
    title('synchronizing TTL post-alignment')
    
elseif ~isempty(ePhysDataPath)&&~isempty(oPhysDataPath) % most likely only ePhys & oPhys data found...
    
    figure(2)
    subplot(211)
    plot(t_E,ttl_E,t_O,ttl_O-1);
    legend('ePhys','oPhys')
    title('synchronizing TTL pre-alignment')
    
    % aligning ePhys & oPhys
    [range_ttlE,range_ttlO]=alignSyncTTL(ttl_E,ttl_O,fs_E);
    ttl_E=ttl_E(range_ttlE);ePhysData=ePhysData(range_ttlE,:);
    ttl_O=ttl_O(range_ttlO,:);oPhysData=oPhysData(range_ttlO,:);
    loco=zeros(size(ttl_O)); %for the final plot

    subplot(212)
    plot(getTime(ttl_E,fs_E),[ttl_E ttl_O-1]);
    legend('ePhys','oPhys')
    title('synchronizing TTL post-alignment')
    
elseif ~isempty(oPhysDataPath)&&~isempty(behaDataPath) % most likely only oPhys & Behavior data found...
    
    figure(2)
    subplot(211)
    plot(t_B,ttl_B+1,t_O,ttl_O-1);
    legend('behavior','oPhys')
    title('synchronizing TTL pre-alignment')
    
    % oPhys & Behavior
    [range_ttlB,range_ttlO]=alignSyncTTL(ttl_B,ttl_O,fs_E);
    ttl_B=ttl_B(range_ttlB);ePhysData=ePhysData(range_ttlB,:);
    ttl_O=ttl_O(range_ttlO,:);oPhysData=oPhysData(range_ttlO,:);
       
    subplot(212)
    plot(getTime(ttl_B,fs_B),[ttl_B ttl_O-1]);
    legend('behavior','oPhys')
    title('synchronizing TTL post-alignment')
    
else
    disp('only one datatype found > no need for data alignment...')
end

%% SAVINNG THE FINAL STRUCTURE
output.ePhysChannels={'lfp1','lfp2','lfp3'};
output.ePhys=ePhysData;
output.oPhysChannels={'gevi','ref','led','nan'};
output.oPhys=oPhysData;
output.locomotion=loco;
output.ttl=ttl_E;
output.samplingRate=fs_E;

figure()
plot(getTime(output.ttl,output.samplingRate),...
    [output.ttl output.locomotion output.ePhys output.oPhys])
legend('ttl','loco','lfp1','lfp2','lfp3','gevi','ref','led','nan')
end

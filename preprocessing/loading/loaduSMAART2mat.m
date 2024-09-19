% analysis file for CELL paper - PV-ASAP3 during rest/run transition

% path='X:\Simon\Calibration\Benchmarking\m0000\20211101\uSmaart_static_000';

% fileInfo=dir(fullfile(path,'**/*mat'));

% load(fullfile(fileInfo.folder,fileInfo.name));
function [data,metadata]=loaduSMAART2mat(filePath, varargin)

% example [data,metadata]=loaduSMAART2mat(filePath,'figure',true)
% metadata.fs=fs;
% metadata.channels=channels;

%% OPTIONS

options.figure=true;
options.verbose=true;

%% UPDATE OPTIONS

if nargin>1
    options=getOptions(options,varargin);
end

%%

if options.verbose
    disp('loading...')
end

dataTemp=load(filePath);

if options.verbose
    disp('done loading!')
end

ttl561=[];ttl561_ts=[];laser561=[];ref=[];xtk=[];
ttl488=[];ttl488_ts=[];laser488=[];vol=[];

fs=double(1/dataTemp.dev3626.demods(1).sample_r_avg{1, 1}.header.gridcoldelta);

nFile=size(dataTemp.dev3626.demods(1).sample_dio_avg,2);

if options.verbose
    nbytes = fprintf('processing %d%% \n', 0);
end

for iFile=1:nFile
    
    % Pull Data from MFLI5442 - 488nm laser line
    ttl488=[ttl488 dataTemp.dev5442.demods(1).sample_dio_avg{1,iFile}.value];
    ttl488_ts=[ttl488_ts dataTemp.dev5442.demods(1).sample_dio_avg{1,iFile}.timestamp];
    
    laser488=[laser488 dataTemp.dev5442.demods(1).sample_r_avg{1,iFile}.value];
    vol=[vol dataTemp.dev5442.demods(2).sample_r_avg{1,iFile}.value];       % demod at 75kHz - blue laser
    % xtkRG=[xtkRG dataTemp.dev5442.demods(3).sample_r_avg{1,iFile}.value]; % demod at 50kHz - green laser
    
    % Pull Data from MFLI3626 - 561nm laser line
    ttl561=[ttl561 dataTemp.dev3626.demods(1).sample_dio_avg{1,iFile}.value];
    ttl561_ts=[ttl561_ts dataTemp.dev3626.demods(1).sample_dio_avg{1,iFile}.timestamp];
    
    laser561=[laser561 dataTemp.dev3626.demods(1).sample_r_avg{1,iFile}.value];
%    SH_20230807 > wrong assigment: red and xtk have been swapped! check for any deletrious mistake in
%     prior analysis - below is how it was prior to correction.
%     ref=[ref dataTemp.dev3626.demods(2).sample_r_avg{1,iFile}.value]; 
%     xtk=[xtk dataTemp.dev3626.demods(3).sample_r_avg{1,iFile}.value];
    xtk=[xtk dataTemp.dev3626.demods(2).sample_r_avg{1,iFile}.value];       % demod at 75kHz - blue laser
    ref=[ref dataTemp.dev3626.demods(3).sample_r_avg{1,iFile}.value];       % demod at 50kHz - green laser
        
    if options.verbose
        fprintf(repmat('\b',1,nbytes))
        nbytes = fprintf('processing %d%% \n', round(iFile/nFile*100));
    end
    
end

%%

if options.figure
    time=getTime(ttl488,fs);
    
    figure;
    subplot(3,2,[1 2])
    plot(time,[ttl488; ttl561]+[0.5; -0.5])
    title('TimeStamps diagnostics')
    
    subplot(323)
    plot(time,ttl488_ts)
    
    subplot(324)
    plot(time,ttl561_ts)
    
    subplot(325)
    histogram(diff(ttl488_ts),100)
    
    subplot(326)
    histogram(diff(ttl561_ts),100)
    
    figure;
    plot(time,vol)
    
end
%% OUTPUT DATA

channels={'laser488',...
    'laser561',...
    'green',...
    'greenXT',...
    'red',...
    'ttl488',...
    'ttl561',...
    'ttl488_ts',...
    'ttl561_ts'};

data(:,1)=laser488;
data(:,2)=laser561;
data(:,3)=vol;
data(:,4)=xtk;
data(:,5)=ref;
data(:,6)=ttl488;
data(:,7)=ttl561;
data(:,8)=ttl488_ts;
data(:,9)=ttl561_ts;

metadata.fs=fs;
metadata.channels=channels;

end

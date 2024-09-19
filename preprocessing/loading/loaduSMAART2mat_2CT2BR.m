% analysis file for CELL paper - PV-ASAP3 during rest/run transition

% path='X:\Simon\Calibration\Benchmarking\m0000\20211101\uSmaart_static_000';

% fileInfo=dir(fullfile(path,'**/*mat'));

% load(fullfile(fileInfo.folder,fileInfo.name));
function [data,metadata]=loaduSMAART2mat_2CT2BR(filePath, varargin)

% example [data,metadata]=loaduSMAART2mat(filePath,'figure',true)
% metadata.fs=fs;
% metadata.channels=channels;

%% OPTIONS

options.figure=true;

%% UPDATE OPTIONS

if nargin>=2
    options=getOptions(options,varargin);
end

%%
disp('loading...')
dataTemp=importdata(filePath);
disp('done loading!')

ttl1=[];ttl1_ts=[];G1=[];xtG1=[];R1=[];
ttl2=[];ttl2_ts=[];G2=[];xtG2=[];R2=[];

fs=double(1/dataTemp.dev3626.demods(1).sample_r_avg{1, 1}.header.gridcoldelta);

nFile=size(dataTemp.dev3626.demods(1).sample_dio_avg,2);

nbytes = fprintf('processing %d%% \n', 0);

for iFile=1:nFile

% Pull Data from MFLI5442 - 488nm laser line
ttl1=[ttl1 dataTemp.dev5442.demods(1).sample_dio_avg{1,iFile}.value];
ttl1_ts=[ttl1_ts dataTemp.dev5442.demods(1).sample_dio_avg{1,iFile}.timestamp];

G1=[G1 dataTemp.dev5442.demods(1).sample_r_avg{1,iFile}.value];
xtG1=[xtG1 dataTemp.dev5442.demods(2).sample_r_avg{1,iFile}.value];
R1=[R1 dataTemp.dev5442.demods(3).sample_r_avg{1,iFile}.value];

% Pull Data from MFLI3626 - 561nm laser line
ttl2=[ttl2 dataTemp.dev3626.demods(1).sample_dio_avg{1,iFile}.value];
ttl2_ts=[ttl2_ts dataTemp.dev3626.demods(1).sample_dio_avg{1,iFile}.timestamp];

G2=[G2 dataTemp.dev3626.demods(1).sample_r_avg{1,iFile}.value];
xtG2=[xtG2 dataTemp.dev3626.demods(2).sample_r_avg{1,iFile}.value];
R2=[R2 dataTemp.dev3626.demods(3).sample_r_avg{1,iFile}.value];

fprintf(repmat('\b',1,nbytes))
nbytes = fprintf('processing %d%% \n', round(iFile/nFile*100));
    
end

%% 

if options.figure
time=getTime(ttl1,fs);

disp(sum(ttl1_ts-ttl2_ts,'all'));

figure;
subplot(3,2,[1 2])
plot(time,[ttl1; ttl2]+[0.5; -0.5])
title('TimeStamps diagnostics')

subplot(323)
plot(time,ttl1_ts)

subplot(324)
plot(time,ttl2_ts)

subplot(325)
histogram(diff(ttl1_ts),100)

subplot(326)
histogram(diff(ttl2_ts),100)

figure;
plot(time,G1)

end
%% OUTPUT DATA

channels={'G1',...
          'xtG1',...
          'R1',...
          'G2',...
          'xtG2',...
          'R2',...
          'ttl1',...
          'ttl2',...
          'ts1',...
          'ts2'};
      
data(:,1)=G1;
data(:,2)=xtG1;
data(:,3)=R1;
data(:,4)=G2;
data(:,5)=xtG2;
data(:,6)=R2;
data(:,7)=ttl1;
data(:,8)=ttl2;
data(:,9)=ttl1_ts;
data(:,10)=ttl2_ts;

metadata.fs=fs;
metadata.channels=channels;

end

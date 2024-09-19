
function [ttl]=findRestRunTransitions(speed,varargin)
% find the mouse state transitino from behavioral speed data (running wheel
% or open field)

%% DEFAULT OPTIONS
options.speedCutOff=2; % in cm/s
options.transCutOff=0.2; % mean speed value to assign 'resting state'
options.maxNumChanges=15;
options.resample=true;
%% UPDATE OPTIONS
if nargin>1
    options=getOptions(options,varargin);
end

%%
rrTransition=speed;
rrTransition(speed>options.speedCutOff)=1;
rrTransition(speed<=options.speedCutOff)=0;

if options.resample
rrTransition=resample(rrTransition,1,100);
end

rrTransition(rrTransition>0.5)=1;
rrTransition(rrTransition<=0.5)=0;

[ipt]=findchangepts(rrTransition,'MaxNumChanges',options.maxNumChanges);

n=length(rrTransition);

% SH-20221201 - changed ipt' to ipt due to dimension error
% SH-20230113 - added row vector check to avoid error
if ~isrow(ipt)
    ipt=ipt';
end

ipt=[1 ipt n]; 
temp=[];
for i=1:numel(ipt)-1
    temp(i)=mean(rrTransition(ipt(i):ipt(i+1)-1));
end

temp(temp<options.transCutOff)=0;temp(temp>=options.transCutOff)=1;

ttl=zeros(1,n);
idx=find(temp>0);
for i=1:numel(idx)
    ttl(ipt(idx(i)):ipt(min(idx(i)+1,numel(ipt))))=1;
end

if options.resample
ttl=resample(ttl,100,1);
end

ttl(ttl>=0.5)=1;
ttl(ttl<0.5)=0;
end

function [output]=getStimEpoch(Data,TTL,Fs,varargin)
% [output]=getStimEpoch(Data,TTL,Fs,'baselinePrePost',2,'getShuffle',true)
% output arrays : time x trials x channels/pixels
%
% OPTIONS:
% options.baselinePrePost=1; % 1sec baseline
% options.getShuffle=false;
% options.preNormalize='median'; % mean %median % none
% options.stimLength=[];
% 
% OUTPUT:
% output.arrayRaw=arrayRaw;
% output.arrayShuffle=arrayShuffle;
% output.stimBand=stimBand;
% output.indexTTLraw=idx_raw;
% output.indexTTLshuffle=idx_shuffle;
% output.options=options;
%
% 20211112 - SH > correct bad increment. Added k=k+1;
% 20211129 - SH > error in shuffle loop > k not reset.
% 20211214 - SH > expand normalization options
% 20221129 - SH > bug in normalization - added omitnan and specify which
% dimension to normalize (frpm median(data(range)) to
% median(data(range,:),'omitnan'));

%% DEFAULT OPTIONS
options.baselinePrePost=1; % 1sec baseline
options.getShuffle=false;
options.preNormalize='median'; % mean %median % none
options.stimLength=[];

%% UPDATE OPTIONS
if nargin>3
    options=getOptions(options,varargin);
end
%%
% assume constant stimulation length
[on,~]=find(diff(TTL) == 1, 1,'first');
[off,~]=find(diff(TTL) == -1, 1,'first');

% find stimulus length
marging=options.baselinePrePost;

if isempty(options.stimLength)
    stimLength=max(round((off-on+1)/Fs,1),1);
else
    stimLength=options.stimLength;
end
if numel(marging)==1
    stimBand=[-abs(marging) marging+stimLength];
    pre=abs(marging); %post=stimBand(2);
elseif numel(marging)==2
    stimBand=marging;
    pre=abs(stimBand(1)); %post=stimBand(2);
else
    error('wrong number of element in baselinePrePost')
end

if istensor(Data)
    options.isTensor=true;
    dim=size(Data);
    Data=reshape(Data,dim(1)*dim(2),dim(3));
    Data=Data';
else
    options.isTensor=false;
end

nChannel=size(Data,2);

% find the 0 before 0-1 transition
[idx_raw,~]=find(diff(TTL) == 1);
if options.getShuffle
    temp=round(1+2*Fs:stimLength*Fs:length(TTL)-2*Fs); %to avoid epoch assignment failure
    temp=temp(randperm(numel(temp)));
    idx_shuffle = sort(temp(1:numel(idx_raw)),'ascend')';
end
% plot([idx_raw idx_shuffle],'o')
nStim=length(idx_raw);

fprintf('%d Cues detected\n',nStim)

%% sort out each epoch for each channels
totalLength=round(diff(stimBand)*Fs);
arrayRaw=zeros(totalLength,nStim,nChannel);
arrayShuffle=zeros(totalLength,nStim,nChannel);

disp('Assigning epoch on Raw data')
k=1;
for i=1:nStim
    try
        start=round(idx_raw(i)-pre*Fs);
        range=(start:start+totalLength-1)';
        switch options.preNormalize
            case 'zscore'
                normRange=(1:floor(pre*Fs)-1)'; % normalized with data prior stimulus
                DataTemp=sh_zscore(Data(range,:),'range',normRange);
            case 'mean'
                normRange=(1:floor(pre*Fs)-1)'; % normalized with data prior stimulus
                DataTemp=Data(range,:)-mean(Data(range(normRange),:),'omitnan');
            case 'median'
                normRange=(1:floor(pre*Fs)-1)'; % normalized with data prior stimulus
                DataTemp=Data(range,:)-median(Data(range(normRange),:),'omitnan');
            case 'none'
                DataTemp=Data(range,:);
        end
        arrayRaw(:,k,:)=DataTemp;
        k=k+1;
    catch
        disp(strcat('fail @ index=',num2str(i)))
    end
end
arrayRaw(:,k:end,:)=[];

if options.getShuffle
    disp('Assigning epoch on Shuffled data')
    k=1;
    for i=1:nStim
        try
            start=round(idx_shuffle(i)-pre*Fs);
            range=(start:start+totalLength-1)';
            switch options.preNormalize
                case 'zcore'
                    normRange=(1:floor(pre*Fs)-1)'; % normalized with data prior stimulus
                    DataTemp=sh_zscore(Data(range,:),'range',normRange);
                case 'mean'
                    normRange=(1:floor(pre*Fs)-1)'; % normalized with data prior stimulus
                    DataTemp=Data(range,:)-mean(Data(range(normRange),:),'omitnan');
                case 'median'
                    normRange=(1:floor(pre*Fs)-1)'; % normalized with data prior stimulus
                    DataTemp=Data(range,:)-median(Data(range(normRange),:),'omitnan');
                case 'none'
                    DataTemp=Data(range,:);
            end
            arrayShuffle(:,k,:)=DataTemp;
            k=k+1;
        catch
            disp(strcat('fail @ index=',num2str(i)))
        end
    end
end
arrayShuffle(:,k:end,:)=[];
%%
if options.isTensor
    trialN=size(arrayRaw,2);
    arrayRaw=permute(arrayRaw,[3 1 2]);
    arrayRaw=reshape(arrayRaw,dim(1),dim(2),[],trialN);
    if options.getShuffle
        arrayShuffle=permute(arrayShuffle(:,1:trialN,:),[3 1 2]);
        arrayShuffle=reshape(arrayShuffle,dim(1),dim(2),[],trialN);        
    end
end

output.arrayRaw=arrayRaw;
output.indexTTLraw=idx_raw;
if options.getShuffle
    output.arrayShuffle=arrayShuffle;
    output.indexTTLshuffle=idx_shuffle;
end
output.stimBand=stimBand;
output.options=options;

end
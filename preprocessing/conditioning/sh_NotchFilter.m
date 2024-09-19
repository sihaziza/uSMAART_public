function [DataF]=sh_NotchFilter(Data,Fsampling,Fcenter,varargin)
% apply a notch filter on a time trace. Data can be a matrix or tensor
% [DataFilt]=sh_NotchFilter(movie,fs,60)
% [DataFilt]=sh_NotchFilter(trace,fs,60,2)

%%
options.deltaFreq=0.5;
options.harmonics=1;
options.order=8;

%% UPDATE OPTIONS
if nargin>3
    options=getOptions(options,varargin);
end

%%
dim=size(Data);

if istensor(Data)
    DataF=reshape(Data,dim(1)*dim(2),dim(3));
elseif dim(1)<dim(2)
    DataF=Data';
else
    DataF=Data;
end

if ~strcmpi(class(DataF),'double')
    DataF=double(DataF);
end

% mind the dimension to which to operate > filtfilt operates on the 1st dim
% > so make sure size(data,1) is time.
% if size(DataF,1)<size(DataF,2)
%     DataF=DataF';
% end
Fsampling=round(Fsampling);

for i=1:options.harmonics
    Fcutoff=double(Fcenter*i);
    if ~isempty(options.deltaFreq)
    lowB=Fcutoff-options.deltaFreq;
    highB=Fcutoff+options.deltaFreq;
    else
    lowB=Fcutoff*0.99;
    highB=Fcutoff*1.01;
    end
    d= designfilt('bandstopiir','FilterOrder',options.order, ...
        'HalfPowerFrequency1',lowB,'HalfPowerFrequency2',highB, ...
        'DesignMethod','butter','SampleRate',Fsampling);
    
    DataF = filtfilt(d,DataF); % mind the dimension along which to operate
end

if istensor(Data)
    DataF=reshape(single(DataF),dim(1),dim(2),dim(3));
elseif dim(1)<dim(2)
    DataF=single(DataF');
else
    DataF=single(DataF);
end

if istensor(Data)
DataF=reshape(DataF,dim(1),dim(2),dim(3));
end
% DataF=cast(DataF,type);
end
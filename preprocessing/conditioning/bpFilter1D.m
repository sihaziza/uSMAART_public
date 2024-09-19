function [dataFilt]=bpFilter1D(data, BW, Fs,varargin)
% Fband=[0.5 100]-band or  [inf 100]-low or [1 ing]-high
% updated SH-20210210
% temporal detrending of fluorescence movie using lowpass method by default. Sampling Rate must be pass through 
% [dataFilt]=bpFilter1D(data, BW, Fs,varargin)
% All options arguments:
% options.methods='lowpass'; %could be exp or spline
% options.samplingRate=[];
% options.lpCutOff=0.5;%in Hz
% options.verbose=1;
% options.plot=true;
% options.diary=false;

%% OPTIONS
options.methods='lowpass'; %could be exp or spline
options.samplingRate=[];
options.lpCutOff=0.5;%in Hz
options.verbose=1;
options.plot=true;
options.diary=false;
options.order=3;

%% UPDATE OPTIONS
if nargin>3
    options=getOptions(options,varargin);
end

%% CHECK OPTIONS VALIDITY
% if strcmp(options.methods,'lowpass') && isempty(options.samplingRate)
%     error('Need to input sampling rate Fs for lowpass method. [data_dtr]=runPhotoBleachingRemoval(data,"samplingRate",Fs)')
% end

%% CORE FUNCTION
d=size(data);
Fs=round(Fs);
filtOrder=options.order;
% check if vector, matrix or tensor
if istensor(data)
temp=reshape(data,d(1)*d(2),d(3)); % often space*space*time
% warning:  filtfilt operates along the first array dimension of x with size greater than 1.
% invert so as to get time on first, row dimension
temp=temp';
elseif d(1)<d(2) % first dimension to be time
    temp=data';
else
    temp=data;
end

if BW(1)==inf
    [b,a]=butter(int32(filtOrder/2),BW(2)/(Fs/2),'low');
elseif BW(2)==inf
    [b,a]=butter(int32(filtOrder/2),BW(1)/(Fs/2),'high');
else
    [b,a]=butter(filtOrder,BW/(Fs/2),'bandpass');
end

% Often single type are passed trough...
temp=double(temp);
vectTrend=filtfilt(b,a,temp);
% plot([temp(290,:)' vectTrend(:,290)+1])
if  numel(d)>2
dataFilt=reshape(vectTrend',d(1),d(2),d(3));
else
dataFilt=vectTrend;
end

dataFilt=cast(dataFilt,'like',data);

% implay(mat2gray(1000.*dataFilt))
end

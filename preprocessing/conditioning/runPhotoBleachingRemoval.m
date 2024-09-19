function [data_dtr]=runPhotoBleachingRemoval(data,varargin)
% temporal detrending of fluorescence movie using lowpass method by default. Sampling Rate must be pass through
% [data_dtr]=runPhotoBleachingRemoval(data,'samplingRate',Fs)
% All options arguments:
% options.methods='lowpass'; %could be exp or spline
% options.samplingRate=[];
% options.lpCutOff=0.1;%in Hz
% options.verbose=1;
% options.plot=true;
% options.diary=false;

%% OPTIONS
options.methods='lowpass'; %could be exp or spline
options.samplingRate=[];
options.lpCutOff=0.1;%in Hz
options.verbose=1;
options.plot=true;
options.diary=false;
options.filterOrder=2;

%% UPDATE OPTIONS
if nargin>=2
    options=getOptions(options,varargin);
end

%% CHECK OPTIONS VALIDITY
if strcmp(options.methods,'lowpass') && isempty(options.samplingRate)
    error('Need to input sampling rate Fs for lowpass method. [data_dtr]=runPhotoBleachingRemoval(data,"samplingRate",Fs)')
end

%% CORE FUNCTION
d=size(data);

% check if vector, matrix or tensor
if numel(d)>2
    temp=reshape(data,d(1)*d(2),d(3));
elseif d(1)>d(2) % need row vector
    temp=data';
else
    temp=data;
end

nanFlag=sum(isnan(temp),'all');

if nanFlag
    disp('NaN values detected... excluding those for subsequent computation')
    idx=isnan(temp);
    temp(idx)=0;
end

Fs=options.samplingRate;
lowF=options.lpCutOff;
order=options.filterOrder;

switch options.methods
    case 'lowpass'
        tic;
        [b,a]=butter(order,lowF/(Fs/2),'low');
        vectTrend=filtfilt(b,a,double(temp'));
        tempCorr=double(temp)./vectTrend'; % normalization to 1
        %         tempCorr=vectTrend'; % no normalization
        toc;
    case 'exp'
        tic;
        tempCorr=zeros(size(temp));
        parfor iPixel=1:d(1)*d(2)
            [~, ~,outputExp] = createFitExp(temp(iPixel,:));
            tempCorr(iPixel,:)=outputExp.residuals;
        end
        toc;
    case 'spline'
        tic;
        tempCorr=zeros(size(temp));
        parfor iPixel=1:d(1)*d(2)
            [~, ~,outputExp] = createFitSpline(temp(iPixel,:));
            tempCorr(iPixel,:)=outputExp.residuals;
        end
        toc;
end

if nanFlag
    disp('NaN values detected... reshaping the matrix back')
    temp(idx)=nan;
end

if  numel(d)>2
    data_dtr=reshape(tempCorr,d(1),d(2),d(3));
else
    data_dtr=tempCorr;
end

data_dtr=single(data_dtr);
end

function [psd,frequency,pow,options]=plotPSD(data,varargin)

% EXAMPLE: [psd,frequency,pow,options]=plotPSD(data,'FrameRate',200,'FreqBand',[0.1 100])
% 'VerboseMessage'  ,true;
% EXAMPLE: [frequency,pow,options]=plotPSD(data,'samplingRate',200,'BW',[0.1 100])
% 'verbose'  ,true;
% 'plotFigure'   ,true;
% 'Savefig'         ,false;
% 'samplingRate'       ,1000;
% 'BW'        ,[0.1 min(options.samplingRate/2,30)];
% 'window'          ,5;
% 'figureHandle'    ,[];
% 'scaleAxis'       ,'linear';
% 'plotAverage'     ,false
% 'overlap'         , 0.9 > between 0:1
% 'nfftFactor'      , 10 

% DEFAULT Options
options.verbose=true;
options.plotFigure=true;
options.saveFig=false;
options.samplingRate=1000;
options.BW=[0.1 min(options.samplingRate/2,30)];
options.window=5;
options.figHandle=[];
options.scaleAxis='linear';%'log2'
options.plotAverage=false;
options.plotColor=[];
options.overlap=0.8;% between 0:1
options.nfftFactor=10;
%%
% USER-DEFINED INPUT OPTIONS
if nargin>1
    options=getOptions(options,varargin);
end
%%
if istensor(data)
    dim=size(data);
    [~,b]=max(dim);
    a=setdiff([1 2 3],b);
    datatemp=reshape(data,dim(b),dim(a(1))*dim(a(2)));
elseif isrow(data)
    datatemp=data';
else
    datatemp=data;
end

Fs=round(options.samplingRate);
win=options.window*Fs;
ovl=round(options.overlap*win);
nfft=10*Fs;

[psd,frequency] =pwelch(datatemp-median(datatemp),win,ovl,nfft,Fs,'onesided');
pow=10*log10(psd);

if options.plotAverage
    pow=mean(pow,2);
end
% k=find(frequency>=band(2),1,'first');
% pow=pow-min(pow(1:k,:),[],1); % to norm all to same noise floor

if options.plotFigure
    if isempty(options.figHandle)
        figure('Name','Power Spectrum Density','DefaultAxesFontSize',16,'color','w');
    else
        figure(options.figHandle)
    end
    switch options.scaleAxis
        case 'linear'
            % plot(log10(frequency),movmin(pow,round(k)/10));
            % band=[0 2];
            % plot(log10(frequency),pow,'linewidth',2)
            trace=pow;%-movmin(pow,20*Nfft/Fs);
            plot(frequency,trace,'linewidth',2);%,'color',options.plotColor)
            xlim(options.BW)
%             title('Power Spectrum Density')
            xlabel('Frequency (Hz)')
            ylabel('Power dB')
%             title('pWelch-estimated Power Spectral Density')
        case 'log10'
            % band=[0 2];
            % plot(log10(frequency),pow,'linewidth',2)
%             options.BW=[-2 log10(Fs/2)];
            trace=pow;%-movmin(pow,20*Nfft/Fs);
            plot(log10(frequency),trace,'linewidth',1.5)
            xlim(log10(options.BW))
            ylabel('Power Spectrum Density')
            xlabel('Frequency (Hz)')
%             title('pWelch-estimated Power Spectral Density')
          case 'log2'
               trace=pow;%-movmin(pow,20*Nfft/Fs);
            plot(log2(frequency),trace,'linewidth',1.5)
            xlim(log2(options.BW))
            ylabel('Power Spectrum Density')
            xlabel('Frequency (Hz)')
%             title('pWelch-estimated Power Spectral Density')
        otherwise
            warning('Scale axis not recognize. Only "linear" or "log" accepted.')
    end
end

if istensor(data)
    pow=reshape(pow,size(pow,1),dim(a(1)),dim(a(2)));
    psd=reshape(psd,size(psd,1),dim(a(1)),dim(a(2)));
end
end
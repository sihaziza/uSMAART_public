function [unmixsource,umxcoeff,options]=umxRLR(source,reference,varargin)
% GOAL: Unmix the two signal using robust linear regression method.
% 
% INPUT:
%     - source: source trace 'sce' to be unmixed from the reference trace 'ref'.
%     - reference:
%     - options (as input struct or variable arguments):
%         * samplingRate: sampling frequency of time trace. (default: 1kHz)
%         * HeartbeatRange: a 2-element-vector. (default: awake > 9-13Hz
%         * Verbose: command window output
% OUTPUT:
%     - unmixsource: unmixed source signal
%     - umxcoeff: unmixing coefficient used to correct signal from reference
%     - options: all options used in the function
% 
% DEPENDENCIES: none
% 
% EXAMPLE 1:
%     [umx,a,options]=umxRLR(V,R); % uses default options
% 
% EXAMPLE 2:
%     [umx,a,options]=umxRLR(V,R,'samplingRate',100);
% 
% CONTACT
% StanfordVoltageGroup@gmail.com
%
% HISTORY
% Created 2020-06-01 by Stanford Voltage Group, Schnitzer Lab
%% CHECK VARIABLES FORMATING
if nargin<2
    cprintf('yellow','This function requires at least 2 vector inputs (source & reference).\n',nargin);
end

if ~isvector(source)||~isvector(reference)
    cprintf('yellow','This function requires 2 vector inputs (source & reference).\n');
end

if mod(nargin-2,2)
    cprintf('red','Variable input arguments work as Name-Value pairs');
end

%% DEFAULT INPUT OPTIONS
options.VerboseMessage=true;
options.VerboseFigure=true;
options.Savefig=false;
options.FigureDir='I:\USERS\Simon\TEST';
options.samplingRate=1000; % default: for TEMPO
options.MouseState=[]; % default: no assumption
options.HeartbeatRange=[];
options.WindowPSD=5;

%% USER-DEFINED INPUT OPTIONS
if nargin>2
    options=getOptions(options,varargin);
end

%% PREPARE OUTPUT OPTIONS
options.FunctionPath=mfilename('fullpath');
options.ExecutionDate=datetime('now');
options.ExecutionDuration=tic;

%% FUNCTION CORE

% % Narrow down the heart beat location
% if strcmpi(options.MouseState,'awake')
%     options.HeartbeatRange=[9 14];
% elseif  strcmpi(options.MouseState,'anesthesia')
%     options.HeartbeatRange=[2 6];
% else
%     options.HeartbeatRange=[]; % default: no assumption
% end

Fs=options.samplingRate;

if isempty(options.HeartbeatRange)
    % Find Hemo peak
    [Fhemo,optsHB]=FindHBpeak(reference,'samplingRate',Fs,'VerboseMessage',false,'VerboseFigure',true);
    Fh=[Fhemo.Location-2 Fhemo.Location+2];
    options.FreqHB=Fhemo.Location;
else
    Fh=options.HeartbeatRange;
    options.FreqHB=Fh;
end

% Data conditioning
M=double([source reference]);
[b,a]=butter(3,Fh/(Fs/2),'bandpass');
Mfilt=filtfilt(b,a,M);

% Run Unmixing
p=robustfit(Mfilt(:,2),Mfilt(:,1));
umxcoeff=p(2);
unmixsource=M(:,1)-umxcoeff.*M(:,2); % compute the residual
[~, MSGID] = lastwarn();
warning('off', MSGID);

%% COMMAND WINDOW OUTPUT

if options.VerboseMessage
    fprintf('Unmixing coefficient is @ %2.3f \n',umxcoeff)
end

M=([source reference unmixsource]);

if options.VerboseFigure
    h=figure('Name','Output Summary for Robust Linear Regression (RLR) method');
    subplot(221)
    plot(M(1:10:end,2),M(1:10:end,1),'+k')
    hold on
    plot(M(1:10:end,2),M(1:10:end,3),'+r')
    hold off
    title('Correlation plot')
    xlabel('Reference')
    
    subplot(222)
    plotPSD([source reference unmixsource],'samplingRate',Fs,...
        'BW',[0.5 50],'figHandle',h);
    title('pWelch PSD plot')
    
    subplot(212)
    t=linspace(0,(length(source)-1)/Fs,length(source))';
    plot(t,M)
    title('Time Trace for all signals')
    xlim([0 2])
    
    if options.Savefig
        if isempty(options.FigureDir)
            disp('Figure not saved. Input save directory.')
        else
            savePDF(h,'RLR Unmixing Summary',options.FigureDir)
        end
    end
end
% options.ExecutionDuration=toc(options.ExecutionDuration);
% options.FindHBpeak=optsHB;
end

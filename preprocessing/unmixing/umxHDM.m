function [unmixsource,umxcoeff,options]=umxHDM(source,reference,varargin)
% GOAL: Minimization of the energy within the heartbeat frequency band. This 
% function works on vector variable only.
% 
% INPUT:
%     - source: input signal to be unmixed 
%     - reference: input reference.
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
%     [umx,a,options]=umxHDM(V,R); % uses default options
% 
% EXAMPLE 2:
%     [umx,a,options]=umxHDM(V,R,'samplingRate',100);
% 
% HISTORY
%     - 2020-06-01 - created by Simon Haziza,PhD (sihaziza@stanford.edu)

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
options.coefRange=[0 4];
%% USER-DEFINED INPUT OPTIONS
if nargin>2
    options=getOptions(options,varargin);
end

%% PREPARE OUTPUT OPTIONS
options.FunctionPath=mfilename('fullpath');
options.ExecutionDate=datetime('now');
options.ExecutionDuration=tic;

%% FUNCTION CORE
Fs=options.samplingRate;

% Narrow down the heart beat location
if isempty(options.HeartbeatRange)
    if strcmpi(options.MouseState,'awake')
        options.HeartbeatRange=[9 14];
    elseif  strcmpi(options.MouseState,'anesthesia')
        options.HeartbeatRange=[2 6];
    else
        options.HeartbeatRange=[2 14]; % default: no assumption
    end
    
    % Find Hemo peak
[Fhemo,optsHB]=FindHBpeak(source,'samplingRate',Fs,'VerboseMessage',true,'VerboseFigure',true);
Fh=[Fhemo.Location-1 Fhemo.Location+1];
options.FreqHB=Fh;

else
    disp('using user defined heartbeat frequency range')
    Fh=options.HeartbeatRange; % a 2-element vector
end
    


% Run optimization
Valpha=(options.coefRange(1):0.01:options.coefRange(2))';
x=source*ones(size(Valpha))';
y=reference*Valpha';
E=bandpower(x-y,Fs,Fh);
[Eval,idx]=min(E);
umxcoeff=Valpha(idx);
unmixsource=source-umxcoeff.*reference;

%% COMMAND WINDOW OUTPUT

if options.VerboseMessage
    fprintf('Unmixing coefficient is @ %2.3f \n',umxcoeff)
end

if options.VerboseFigure
    h=figure('Name','Output Summary for HemoDynamic Minimization (HDM) method');
    subplot(221)
    plot(Valpha,E','k','linewidth',2)
    hold on
    plot([umxcoeff umxcoeff],[Eval Eval],'+r','linewidth',2,'markersize',10)
    hold off
    title('Evolution of Hemodynamic peak power')
    xlabel('unmixing coefficient')
    
    M=[reference source unmixsource];
    M=M-mean(M);
    subplot(222)
    plotPSD(M,'samplingRate',Fs,...
        'BW',options.HeartbeatRange,'figHandle',h);
    title('pWelch PSD plot')
    
    subplot(212)
    t=linspace(0,(length(source)-1)/Fs,length(source))';
    plot(t,M)
    title('Time Trace for Reference and Source signals')
    xlim([0 10*1/Fh(1)])
    
    if options.Savefig
        if isempty(options.FigureDir)
            disp('Figure not saved. Input save directory.')
        else
            savePDF(h,'HDM Unmixing Summary',options.FigureDir)
        end
    end 
end
%     options.ExecutionDuration=toc(options.ExecutionDuration);
%     options.FindHBpeak=optsHB;
end

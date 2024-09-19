function [umxSource,umxCoeff,options]=unmixing1D(source,reference,varargin)
% GOAL: Unmix the two signal using various method.
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
% options.message=true;
% options.figure=true;
% options.Savefig=false;
% options.FigureDir='I:\USERS\Simon\TEST';
% options.samplingRate=1000; % default: for TEMPO
% options.MouseState=[]; % default: no assumption
% options.HeartbeatRange=[];
% options.WindowPSD=5;
% options.method='HDM';
% options.freqHemo=[];
%
% DEPENDENCIES: if ICA method is used > https://research.ics.aalto.fi/ica/fastica/
%
% EXAMPLE 1:
%     [umx,a,options]=unmixing1D(V,R); % uses default options
%
% EXAMPLE 2:
%     [umx,a,options]=unmixing1D(V,R,'method','pca');
%
% CONTACT
% StanfordVoltageGroup@gmail.com
%
% HISTORY
% Created 2020-06-01 by Stanford Voltage Group, Schnitzer Lab


%% CHECK VARIABLES FORMATING
if nargin<2
    error('This function requires at least 2 vector inputs (source & reference).\n');
end

if ~isvector(source)||~isvector(reference)
    error('This function requires 2 vector inputs (source & reference).\n');
end

if mod(nargin-2,2)
    error('Variable input arguments work as Name-Value pairs.\n');
end

%% DEFAULT INPUT OPTIONS
options.message=true;
options.figure=true;
options.Savefig=false;
options.FigureDir='I:\USERS\Simon\TEST';
options.samplingRate=1000; % default: for TEMPO
options.MouseState=[]; % default: no assumption
options.HeartbeatRange=[];
options.WindowPSD=5;
options.method='HDM';
options.freqHemo=[];
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

fs=options.samplingRate;
mState=options.MouseState;
vMess=options.message;
vFig=options.figure;
Fhemo=options.freqHemo;
HeartbeatRange=options.HeartbeatRange;
coefRange=options.coefRange;
% Run Unmixing
switch lower(options.method)
    case 'cnv'
        if options.message
            cprintf('yellow','Unmixing with method [%s].\n','HDM');
        end
        [umxSource,umxCoeff]=umxCONV(source,reference,fs...
            ,'MouseState',mState,'HeartbeatRange',HeartbeatRange,...
            'coefRange',coefRange,'VerboseMessage',vMess,'VerboseFigure',vFig);
    case 'hdm'
        if options.message
            cprintf('yellow','Unmixing with method [%s].\n','HDM');
        end
        [umxSource,umxCoeff,umxOptions]=umxHDM(source,reference,...
            'samplingRate',fs,'MouseState',mState,'HeartbeatRange',HeartbeatRange,...
            'coefRange',coefRange,'VerboseMessage',vMess,'VerboseFigure',vFig);
    case 'rlr'
        if options.message
            cprintf('yellow','Unmixing with method [%s].\n','RLR');
        end
        [umxSource,umxCoeff,umxOptions]=umxRLR(source,reference,...
            'samplingRate',fs,'MouseState',mState,'HeartbeatRange',HeartbeatRange,...
            'VerboseMessage',vMess,'VerboseFigure',vFig);
    case 'pca'
        if options.message
            cprintf('yellow','Unmixing with method [%s].\n','PCA');
        end
        [umxSource,umxCoeff,umxOptions]=umxPCA(source,reference,...
            'samplingRate',fs,'freqHemo',Fhemo,'MouseState',mState,...
            'VerboseMessage',vMess,'VerboseFigure',vFig);
    case 'ica'
        if options.message
            cprintf('yellow','Unmixing with method [%s].\n','ICA');
        end
        [umxSource,umxCoeff,umxOptions]=umxICA(source,reference,...
            'samplingRate',fs,'MouseState',mState,...
            'VerboseMessage',vMess,'VerboseFigure',vFig);
end

umxSource=umxSource-mean(umxSource);
options.umxOptions=umxOptions;
options.ExecutionDuration=toc(tic);
% umxSource=source-bpFilter1D(reference,[inf 30],fs)
end

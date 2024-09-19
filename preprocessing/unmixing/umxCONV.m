function [umx,ref_clean,w1,options]=umxCONV(source, reference, fs, varargin)

% convolutional unmixing (adapted from Vasily method).
% [umx,ref_clean,w1,options]=umxCONV(source, reference, fs, varargin)

%% CHECK VARIABLES FORMATING
if nargin<3
    cprintf('yellow','This function requires at least 2 vector inputs (source & reference) and the sampling frequency.\n',nargin);
end

if ~isvector(source)||~isvector(reference)
    cprintf('yellow','This function requires 2 vector inputs (source & reference).\n');
end

%% DEFAULT INPUT OPTIONS

% options for filter estimation:
options.epoch=1; %time segment in sec
options.overlap=0.9; % fraction
options.df_reg=5; 
options.eps=1e-8;
options.max_phase=pi; 
options.VerboseMessage=true;
options.figure=false;
options.Savefig=false;
options.MouseState=[]; % default: no assumption
options.HeartbeatRange=[];
options.WindowPSD=5;
options.coefRange=[0 1];
options.mirrorPadding=true;
options.fref=[];
options.flim_max=[];

%% USER-DEFINED INPUT OPTIONS
if nargin>3
    options=getOptions(options,varargin);
end

%% highpass isn't good imho

dn = ceil(fs*options.epoch);
dn_overlap = round(dn*(1-options.overlap));
df_n=ceil(options.df_reg*options.epoch);
reg_func=@(z) smoothdata(z, 'movmedian', df_n);

source=source-mean(source);
reference=reference-mean(reference);

[w1, ~] = estimateFilter(source, reference, dn, dn_overlap, options.eps, reg_func, options.max_phase);
% [w1, ~] = estimateFilterReg(source, reference, dn, dn_overlap,'max_phase', options.max_phase,'fref',options.fref,'flim_max',options.flim_max);

if options.mirrorPadding
    ref_clean=[(reference(1:dn)); reference; (reference(end-dn+1:end))];
else 
    ref_clean=reference;
end

ref_clean = conv(ref_clean, w1, 'same');
ref_clean = [ref_clean(2:end); ref_clean(end)];

if options.mirrorPadding
    % trim back to the original size
    ref_clean=ref_clean(dn+1:end-dn,:);
end

umx=source-ref_clean;

%% plot summary figure

M=([source ref_clean umx]);
M=M-mean(M);

if options.figure
    h=figure('Name','Output Summary for Robust Linear Regression (RLR) method');
    subplot(221)
    plot(M(1:10:end,2),M(1:10:end,1),'+k')
    hold on
    plot(M(1:10:end,2),M(1:10:end,3),'+r')
    hold off
    title('Correlation plot')
    xlabel('Reference')
    
    subplot(222)
    plotPSD(M,'samplingRate',fs,...
        'BW',[0.5 50],'figHandle',h);
    title('pWelch PSD plot')
    
    subplot(212)
    plot(getTime(M,fs),M)
    title('Time Trace for all signals')
%     xlim([0 10])
    
    if options.Savefig
        if isempty(options.FigureDir)
            disp('Figure not saved. Input save directory.')
        else
            savePDF(h,'RLR Unmixing Summary',options.FigureDir)
        end
    end
end


end


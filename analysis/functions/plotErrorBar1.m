% ----------------------------------------------------------------------- %
% Function plot_areaerrorbar plots the mean and standard deviation of a   %
% set of data filling the space between the positive and negative mean    %
% error using a semi-transparent background, completely customizable.     %
%                                                                         %
%   Input parameters:                                                     %
%       - data:     Data matrix, with rows corresponding to observations  %
%                   and columns to samples.                               %
%       - options:  (Optional) Struct that contains the customized params.%
%           * options.figHandle:       Figure handle to plot the result.     %
%           * options.color_area:   RGB color of the filled area.         %
%           * options.color_line:   RGB color of the mean line.           %
%           * options.alpha:        Alpha value for transparency.         %
%           * options.line_width:   Mean line width.                      %
%           * options.x_axis:       X time vector.                        %
%           * options.error:        Type of error to plot (+/-).          %
%                   if 'std',       one standard deviation;               %
%                   if 'sem',       standard error mean;                  %
%                   if 'var',       one variance;                         %
%                   if 'c95',       95% confidence interval.              %
% ----------------------------------------------------------------------- %
%   Example of use:                                                       %
%       data = repmat(sin(1:0.01:2*pi),100,1);                            %
%       data = data + randn(size(data));                                  %
%       plot_areaerrorbar(data);                                          %
% ----------------------------------------------------------------------- %
%
%   adapted from:  Victor Martinez-Cagigal                                      %
%   Date:    30/04/2018                                                   %
%   E-mail:  vicmarcag (at) gmail (dot) com                               %
% ----------------------------------------------------------------------- %
function plotErrorBar1(data, varargin)
% EXAMPLE: [frequency,pow,options]=plotPSD(data,'FrameRate',200,'FreqBand',[0.1 100])
% 'VerboseMessage'  ,true;
% 'plotFigure'   ,true;
% 'Savefig'         ,false;
% 'FrameRate'       ,1000;
% 'FreqBand'        ,[0.1 min(options.FrameRate/2,30)];
% 'Window'          ,5;
% 'figureHandle'    ,[];
% 'scaleAxis'       ,'linear';
% 'plotAverage'     ,false
% 'overlap'         , 0.9 > between 0:1
% 'nfftFactor'      , 10

% DEFAULT Options
options.figHandle     = [];
options.color_area = [0.5 0.5 0.5];    % dark theme
options.color_line = [0 0 0];
options.alpha      = 0.5;
options.line_width = 2;
options.line_style='-';
options.error      = 'sem';
options.x_axis=[];
%%
% USER-DEFINED INPUT OPTIONS
if nargin>1
    options=getOptions(options,varargin);
end
%%
% assume samples in column
if isrow(data)
    data=data;
else
    data=data';
end

if isempty(options.x_axis)
    options.x_axis = 1:size(data,2);
else
    if isrow(options.x_axis)
        options.x_axis=options.x_axis;
    else
        options.x_axis=options.x_axis';
    end
end

% Computing the mean and standard deviation of the data matrix
data_mean = mean(data,1);
data_std  = std(data,0,1);

% Type of error plot
switch(options.error)
    case 'std', error = data_std;
    case 'sem', error = (data_std./sqrt(size(data,1)));
    case 'var', error = (data_std.^2);
    case 'c95', error = (data_std./sqrt(size(data,1))).*1.96;
    case 'c99', error = (data_std./sqrt(size(data,1))).*2.576;
end

% Plotting the result
if ~isempty(options.figHandle)
    figure(options.figHandle)
else
    figH=figure('DefaultAxesFontSize',12,'color','w');
end
x_vector = [options.x_axis, fliplr(options.x_axis)];
y_vector=[data_mean+error,fliplr(data_mean-error)];
% remove inf value (e.g. log10 PSD which give inf at 0Hz
TFinf = isinf(y_vector)+isinf(x_vector);TFnan = isnan(y_vector)+isnan(x_vector);
TF=TFnan+TFinf;TF=~TF;
% if sum(TF,'all')>0
x_vector=x_vector(TF);y_vector=y_vector(TF);
% end
patch = fill(x_vector,y_vector, options.color_area);
set(patch, 'edgecolor', 'none');
set(patch, 'FaceAlpha', options.alpha);
hold on;
plot(options.x_axis, data_mean,'linestyle',options.line_style, 'color', options.color_line, ...
    'LineWidth', options.line_width);
hold off;

end
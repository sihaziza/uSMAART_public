function [speed]=getMouseSpeedDLC(filename,varargin)
% estimate speed on mouse tracking data from DeeplabCut.
% usually only tracking head (entree 1:3) and tailbase (entree 4:6).

%% DEFAULT OPTIONS
options.bodypart='tail'; % 'head'
options.pixSize=0.114; % in cm; 
options.dispCutOff=10; % in pixel
options.frameRate=20; % standard FPS for behavioral video
options.win=1; % in sec - smoothing window for speed

%% UPDATE OPTIONS
if nargin>1
    options=getOptions(options,varargin);
end
%%

beh=h5read(filename,'/df_with_missing/table/');
bFs=20;
pix=0.114; %cm
pix=options.pixSize; % in cm; 
bFs=options.frameRate; 

head_x=beh.values_block_0(1,:);
head_y=beh.values_block_0(2,:);
tail_x=beh.values_block_0(4,:);
tail_y=beh.values_block_0(5,:);

n=length(head_x);

switch options.bodypart
    case 'tail'
        tx1=tail_x(1:end-1);
        tx2=tail_x(2:end);
        ty1=tail_y(1:end-1);
        ty2=tail_y(2:end);
    case 'head'
        tx1=head_x(1:end-1);
        tx2=head_x(2:end);
        ty1=head_y(1:end-1);
        ty2=head_y(2:end);
end

D_pre=((ty2-ty1).^2+(tx2-tx1).^2).^0.5;
T=getTime(head_x,bFs);
T_corr=getTime(D_pre,bFs);
T_corr(D_pre>options.dispCutOff)=[];
D=D_pre;D(D>options.dispCutOff)=[];

% interpolate outlier datapoints.
D_int=interp1(T_corr,D,T,'pchip');

% estimate speed
speed=D_int*pix./(1/bFs);

% smooth speed 
win=options.win*bFs;
speed=movmedian(speed,win);
end
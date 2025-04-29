function [spatialFootprint, timeTrace,boundbox]=getUnitsROI(h5Path,Fs,varargin)
% to extract time traces based on ROI without loading the full movie
% [units, trace]=getUnitsROI(h5Path,Fs)
% exemple 1:
% [~, trace, boundbox]=getUnitsROI(h5Path,Fs);
% exemple 2:
% [~, trace]=getUnitsROI(h5Path,Fs,'boundbox',boundbox);

% Step1: load 1st second and generate the average frame to crop from
% Step2: get the average pixel trace from the h5 file
% drawcircle
% imcrop

%% DEFAULT INPUT OPTIONS
options.boundbox=[];
options.figure=true;

%% USER-DEFINED INPUT OPTIONS
if nargin>2
    options=getOptions(options,varargin);
end

%% LOAD RAW DCMIG DATA
% disp('Starting Neuron Demixing (PCAICA)')

[~,~,ext]=fileparts(h5Path);

if strcmpi(ext,'.h5')
    disp('h5 file detected')
else
    error('not a h5 or Tiff file')
end

meta=h5info(h5Path);
if size(meta.Datasets)==1 % in case only 'mov' present - SH 20220321
    dim=meta.Datasets(1).Dataspace.Size;
else % assume dataset 'mov' is second, 'fps' is first - SH 20220105
    dim=meta.Datasets(2).Dataspace.Size;
end
mx=dim(1);my=dim(2);numFrame=dim(3);
% dataset=strcat(meta.Name,meta.Datasets.Name);
dataset='/mov';

movie=h5read(h5Path,dataset,[1 1 1],[mx my Fs]);
MEAN_PROJECTION=mean(movie,3);
f0=bpFilter2D(MEAN_PROJECTION,25,1,'parallel',false);

if isempty(options.boundbox)
    crop = get_circular_mask(f0);
    % get_rectangular_mask : to do SH-20220321
    image=f0.*crop;
    [~, boundbox] = autoCropImage(image);
else
    boundbox=options.boundbox;
end

spatialFootprint=h5read(h5Path,dataset,[boundbox(2) boundbox(1) 1],[boundbox(4) boundbox(3)  numFrame]);

if options.figure
figure;
subplot(121)
imshow(MEAN_PROJECTION,[])
subplot(122)
imshow(mean(spatialFootprint,3),[])
end

timeTrace=getPointProjection(spatialFootprint);

% figure(10)
% time=getTime(timeTrace,Fs);
% plot(time,trace)


% DESCRIPTION
%
% SYNTAX
%
% EXAMPLE
%
% CONTACT: StanfordVoltageGroup@gmail.com

%% GET DEFAULT OPTIONS

% options.verbose=true;
% options.boundbox=[];
% options.threshold=0.1; % remove 10% of above background level
% outputdir=fullfile(folderPath,'DemixingPCAICA',fname);
% if ~exist(outputdir,'dir')
%     mkdir(outputdir)
% end
end
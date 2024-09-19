function preDLCconditioning(path)
% this script is batch processing .avi behavioral videos before applying
% DLC to track tail & head body parts. input argument allPath can be just
% one path or a cell array of data path for batch processing.
% the first step is 'cropping' and will ask for user input.

nFile=numel(path);

if nFile==1
    allPath{1}=path;
else
    allPath=path;
end

for i=1:nFile
    file=dir(fullfile(allPath{i},'*/*.avi'));
    
    filePath=fullfile(file.folder,file.name);
    v = VideoReader(filePath);
    
    disp('User should draw the ROI')
    image=readFrame(v);
    
    disp('converting to grayscale...')
    image=rgb2gray(image);
    
    [~,rect(i,:)] = imcrop(image);
    rect=round(rect,0);
    
    % change width and height to the next following multiple of 2,
    % i.e. odd numbers, for MP4 saving (H.264 codec).
    corr=mod(rect,2); rect=rect-abs(corr-1);
end

close all

parfor i=1:nFile    
    file=dir(fullfile(allPath{i},'*/*.avi'));
    
    filePath=fullfile(file.folder,file.name);
    
    [video]=loadAVI(filePath,'crop',true,'cropRect',rect(i,:),'binning',2);
    
    disp('saving movie as mp4')
    [~,name,~]=fileparts(filePath);tic;
    renderMovie(video,fullfile(allPath{i},[name '_trim_bin']),20,'MPEG-4');toc;    
end

end


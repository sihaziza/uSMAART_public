function [outputDLC]=estimateBehavioralState(path)

file=dir(fullfile(path,'*.h5'));
fileMP4=dir(fullfile(path,'*.mp4'));

fileName=fullfile(file.folder,file.name);
disp(fileName);

% [folder,name]=fileparts(fileName);
% movieFile=fullfile(folder,[name '_labeled.mp4']);
movieFile=fullfile(fileMP4.folder,fileMP4.name);
pixDimension=getPixelCalibration(movieFile);

[speed_head]=getMouseSpeedDLC(fileName,'bodypart','head','pixSize',pixDimension);
[speed_tail]=getMouseSpeedDLC(fileName,'bodypart','tail','pixSize',pixDimension);

outputDLC.speed_head=speed_head;
outputDLC.speed_tail=speed_tail;
outputDLC.fs=20;

T=getTime(speed_head,outputDLC.fs);
plot(T,[speed_head; speed_tail])
end
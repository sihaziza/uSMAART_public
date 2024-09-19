function pixDimension=getPixelCalibration(movieFile)
% calibrate pixel dimension

disp('calibrating pixel size...')

% path='Y:\Simon\fiber_1BR1CT\OpenField\m404\20231122\meas00';
% fileMP4=dir(fullfile(path,'*.mp4'));
% movieFile=fullfile(fileMP4.folder,fileMP4.name);

v = VideoReader(movieFile);
numFrame=round(v.Duration*v.FrameRate,0);
video=zeros(v.Height,v.Width,numFrame,'uint8');

tic;
k=1;
p=1;
f = waitbar(0,'Please wait... loading data');
while hasFrame(v)
    frame = readFrame(v);
    video(:,:,k)=rgb2gray(frame);
    k=k+1;
    if k>p*(numFrame-1)/10
        waitbar(0.1*p,f,strcat(num2str(10*p),'% loaded'));
        p=p+1;
    end
    if k>100
        break
    end
end
close(f)
toc;


frame=mean(video,3);
yproj=squeeze(mean(video,2));%imshow(yproj,[])
yproj=rescale(yproj,0,1);
xproj=squeeze(mean(video,1));%imshow(flipud(xproj'),[])
xproj=rescale(xproj,0,1);

% figure
% plot(mean(xproj,2))
% hold on
% plot(mean(yproj,2))
% hold off

thres=0.25;
temp=rescale(mean(xproj,2),0,1);
length=find(temp>thres,1,'last')-find(temp>thres,1,'first');
temp=rescale(mean(yproj,2),0,1);
height=find(temp>thres,1,'last')-find(temp>thres,1,'first');
% manually calibrating thres using video0001 11-26-36_trim_binDLC_resnet50_openFieldAug31shuffle1_500000_filtered_labeled.mp4'
% > 263 & 138

pixDimension=mean([30/length 18/height]);

end
function [cleanID]=cleanSWR(allRipples)


% figure position for home workstation [1383,953,1920,970] / lab workstation [6,1241,1920,970]
% surfacePro [77,148,1262,676]
options.skipCleanFile=true;
options.positionFig=[6,1241,1920,970];%1383,953,1920,970];%[6,1241,1920,970];%%[1005,969,1920,970];1005,969,1920,963
options.frameRate=2000;
options.filePath=[];
options.waitForUser=true;

%% UPDATE OPTIONS
if nargin>1
    options=getOptions(options,varargin);
end

%%

nSWR=min(size(allRipples));

allRipples=bpFilter1D(allRipples,[25 inf],options.frameRate);
fs=options.frameRate;
time=getTime(allRipples(:,1),fs);

step=1;
p=0;
garbageID=[];
while p<=size(allRipples,2)
    figure('name','Press any key to enter your favorite cell number','defaultaxesfontsize',12,'color','w','Position',options.positionFig);
    nUnits=min(25,size(allRipples,2)-p);
    for i=1:nUnits
        subplot(5,5,i)
        plot(time,allRipples(:,p+i))
        title(num2str(p+i))      
        axis tight
        axis off
    end
    xlabel('Time (s)')
    if options.waitForUser
        prompt = {'which ripples to garbage? [as a list, just space e.g. 1 5 8 9 15'};
        dlgtitle = 'Input';
        dims = [2 50];
        pause 
        answer = inputdlg(prompt, dlgtitle,dims);
        garbageID=[garbageID str2num(answer{1})];
    end
    step=step+1;
    p=25*(step-1);
end

cleanID=setdiff(1:nSWR,garbageID);

close all

% % double check with the user
% figH=figure('defaultaxesfontsize',12,'color','w','Name','summary of selected cells','Position',options.positionFig);
% nUnits=max(10,numel(cellList));
% for i=1:numel(cellList)
%    
%     subplot(nUnits,10,[10*(i-1)+3 10*i])
%     plot(time,100*temporal(:,cellList(i)))
%     xlim([0 time(end)])
%     ylabel('dF/F (%)')
%     axis tight
%     title('Press any key to enter your favorite cell number')
% end
% xlabel('Time (s)')
% 
% answer = questdlg('Are you happy with you cell choice? (press any key to continue)', ...
%     'Cell selection Summary', ...
%     'Yes Happy','No do it again','Yes Happy');
% 
% % Handle response
% switch answer
%     case 'Yes Happy'
%         disp('outputting results..')
% %         close all
%     case 'No do it again'
%         close all
%         disp('Lets run cell-check again')
%         [cellList,figH]=extractCheckCellManual(output);
% end


end
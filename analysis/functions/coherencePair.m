function [C,P,F]=coherencePair(wave1, wave2,Fs,varargin)
% [C,P,F]=coherencePair(wave1, wave2,Fs,varargin)
% wave1 and wave2 is always single vector
% options.plot=true;
% options.window=2;
% options.channelName=[];
% options.figHandle=[];
% options.BW=[0.5 50];

options.plot=true;
options.window=2;
options.channelName=[];
options.figHandle=[];
options.BW=[0.5 50];
options.logAxis=false;
options.plotPhase=false;

%% UPDATE OPTIONS
if nargin>3
    options=getOptions(options,varargin);
end

%%

wave1=wave1-median(wave1);
wave2=wave2-median(wave2);

nW1=min(size(wave1));
nW2=min(size(wave2));

% global Fs
window = options.window*Fs;
overlap=round(0.8*options.window*Fs);
Nfft=10*Fs;

C=[];P=[];
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
for i=1:nW1
    %     for j=min(max(nW2-1,1),i+1):max(nW2-1,1) % if no diagonal wanted
    for j=min(max(nW2-1,1),i):max(nW2-1,1)
        if options.plotPhase
        [P(:,i,j)]= cpsd(wave1(:,i),wave2(:,j),window,overlap,Nfft,Fs);
        end
        [C(:,i,j),F] = mscohere(wave1(:,i),wave2(:,j),window,overlap,Nfft,Fs);
    end
end

if options.plot
    % Plotting the result
    if ~isempty(options.figHandle)
        figure(options.figHandle)
    else
        figH=figure('DefaultAxesFontSize',12,'color','w');
    end
    for i=1:nW1
        for j=min(max(nW2-1,1),i+1):max(nW2-1,1)
                        subplot(nW1,max(nW2-1,1),j+(i-1)*(max(nW2-1,1)-1))
                        if options.plotPhase
                        yyaxis left
                        plot(F,180.*angle(squeeze(P(:,i,j)))./pi,'linewidth',2)
                        xlabel('Hz')
                        ylabel('\Theta(f)')
                        title('Cross Spectrum Phase')
                        xlim(options.BW)
                        ylim([-180 180])
            
                        yyaxis right
                        end
            if options.logAxis
                plot(F,10*log10(squeeze(C(:,i,j))),'linewidth',2)
                
            else
                plot(F,squeeze(C(:,i,j)),'linewidth',2)
            end
            title('Magnitude-Squared Coherence')
            ylabel('Coherence')
            xlabel('Hz')
            xlim(options.BW)
            %                 ylim([0 1])
            if ~isempty(options.channelName)
                title([options.channelName{i} options.channelName{nW1+j}])
            end
        end
    end
end

end


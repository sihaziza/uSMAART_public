function BW_output=convertBWtoF(BW_input,F)
% BW_input and BW_output are 2 element vector
% F is the frequency vector, assumed to be in increasing order
% BW_output=convertBWtoF(BW_input,F)

% check if F is sorted in ascending order.
sortingOrder=sign(F(end)-F(1));
if  sortingOrder<0
    [BW_output(1),~,~]=find(F<=BW_input(1),1,'first');
    [BW_output(2),~,~]=find(F>=BW_input(2),1,'last');
else
    [BW_output(1),~,~]=find(F>=BW_input(1),1,'first');
    [BW_output(2),~,~]=find(F<=BW_input(2),1,'last');
end

BW_output=sort(BW_output);

end
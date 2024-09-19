function [rgb]=geviColor(name)
% assign color code to various GEVI
% ace pace vnm pacer ref

% Ace: A8D04E (RGB: 168, 208, 78)
% PACE: 8AD4E4 (RGB: 138, 212, 228)
% VARNAM: F58F8F (RGB: 245, 143, 143)
% PACER: EEA1C6 (238, 161, 198)

% round #2
%pace: 52A0A8
% ace: 9AAB3A
% round #3
% Ace: 9AAB3A
% pAce: 53A0A8
% VARNAM/pAcer: B04F4F F59090 HEX2RGB('B14F50') 
% ref: F3B083

switch name
    case 'ace'
        rgb=[154, 171, 58]./255;
    case 'pace'
        rgb=[ 83, 160, 168]./255;
    case 'vnm'
        rgb=[177    79    80]./255;
    case 'pacer'
        rgb=[245   144   144]./255;
     case 'ref'
        rgb=[243 176 131]./255;
end

end
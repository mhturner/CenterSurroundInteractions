clear all; close all; clc;

centerColor = [127, 191, 123] ./ 255;
surroundColor = [175, 141, 195] ./ 255;

subunitSigma = 12;
subunitSurroundSigma = 150;
surroundStrength = 0.1;

centerSigma = 50;
subunitSpacing = 20;


x = -400:400;

centerProfile = exp((-(x).^2)./(2*centerSigma^2));
centerProfile = centerProfile./(max(centerProfile)); %peaks at 1

surroundProfile = exp((-(x).^2)./(2*subunitSurroundSigma^2));
surroundProfile = -surroundStrength * surroundProfile./(max(surroundProfile));

figure; clf; fig1=gca;
set(fig1,'XScale','linear','YScale','linear')
set(0, 'DefaultAxesFontSize', 12)
set(get(fig1,'XLabel'),'String','Pos')
set(get(fig1,'YLabel'),'String','Res')
set(gcf, 'WindowStyle', 'docked')

offsets = -400:subunitSpacing:400;
for pp = 1:length(offsets)
    if pp == 21
        color = centerColor;
    else
        color = [0.5 0.5 0.5];
    end
offset = offsets(pp);
ind = find(offset == x);
ampMult = centerProfile(ind);

amp = exp((-(offset-x).^2)./(2*subunitSigma^2));
amp = ampMult .* (amp./(max(amp))); %integrates to 1, then mult up

amp_s = exp((-(offset-x).^2)./(2*subunitSurroundSigma^2));
amp_s = surroundStrength*ampMult .* (amp_s./(max(amp_s))); %integrates to 1, then mult up

addLineToAxis(x,amp-amp_s,['sub',num2str(pp)],fig1,color,'-','none')

end
addLineToAxis(x,centerProfile,'center',fig1,'k','-','none')




figure; clf; fig2=gca;
set(fig2,'XScale','linear','YScale','linear')
set(0, 'DefaultAxesFontSize', 12)
set(get(fig2,'XLabel'),'String','Pos')
set(get(fig2,'YLabel'),'String','Res')
set(gcf, 'WindowStyle', 'docked')

offsets = -400:subunitSpacing:400;
for pp = 1:length(offsets)
    if pp == 21
        color = centerColor;
    else
        color = [0.5 0.5 0.5];
    end
offset = offsets(pp);
ind = find(offset == x);
ampMult = centerProfile(ind);

amp = exp((-(offset-x).^2)./(2*subunitSigma^2));
amp = ampMult .* (amp./(max(amp))); %integrates to 1, then mult up


addLineToAxis(x,amp,['sub',num2str(pp)],fig2,color,'-','none')


end
addLineToAxis(x,centerProfile,'center',fig2,'k','-','none')
addLineToAxis(x,surroundProfile,'surround',fig2,surroundColor,'-','none')

makeAxisStruct(fig1,'RFtoon_subSurround' ,'RFSurroundFigs')
makeAxisStruct(fig2,'RFtoon_subClassic' ,'RFSurroundFigs')



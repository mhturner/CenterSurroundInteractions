function doNatImageAdditivityAnalysis(node,varargin)
    ip = inputParser;
    expectedMetrics = {'integrated','peak'};
    ip.addRequired('node',@(x)isa(x,'edu.washington.rieke.jauimodel.AuiEpochTree'));
    addParameter(ip,'metric','integrated',...
        @(x) any(validatestring(x,expectedMetrics)));
    addParameter(ip,'figureID',[],@ischar);

    ip.parse(node,varargin{:});
    node = ip.Results.node;
    
    metric = ip.Results.metric;
    figureID = ip.Results.figureID;

    figure; clf;
    fig1=gca; %example C trace
    set(fig1,'XScale','linear','YScale','linear')
    set(0, 'DefaultAxesFontSize', 12)
    set(get(fig1,'XLabel'),'String','Time (s)')

    figure; clf;
    fig2=gca; %example S trace
    set(fig2,'XScale','linear','YScale','linear')
    set(0, 'DefaultAxesFontSize', 12)
    set(get(fig2,'XLabel'),'String','Time (s)')
    
    figure; clf;
    fig3=gca; %example C + S trace and linear sum trace
    set(fig3,'XScale','linear','YScale','linear')
    set(0, 'DefaultAxesFontSize', 12)
    set(get(fig3,'XLabel'),'String','Time (s)')
    
    figure; clf;
    fig4=gca; %summary plot: measured vs linear sum
    set(fig4,'XScale','linear','YScale','linear')
    set(0, 'DefaultAxesFontSize', 12)
    set(get(fig4,'XLabel'),'String','R(C) + R(S)')
    set(get(fig4,'YLabel'),'String','R(C + S)')
    
    cellInfo = getCellInfoFromEpochList(node.epochList);
    recType = getRecordingTypeFromEpochList(node.epochList);

    colors = pmkmp(4);
    for ii = 1:node.children.length %loop over images
        imageNode = node.children(ii);
        responses = zeros(1,4); %c, s, cs, c + s
        ct = 0;
        for ll = 1:imageNode.children.length %loop over patch locations
            locationNode = imageNode.children(ll);
            if locationNode.children.length < 3
                continue
            end
            ct = ct + 1;
            centerResponse = getMeanResponseTrace(locationNode.childBySplitValue('Center').epochList,recType);
            surroundResponse = getMeanResponseTrace(locationNode.childBySplitValue('Surround').epochList,recType);
            centerSurroundResponse = getMeanResponseTrace(locationNode.childBySplitValue('Center-Surround').epochList,recType);
            linearSumResponse = centerResponse.mean + surroundResponse.mean;

            tempCenter = getResponseAmplitudeStats(locationNode.childBySplitValue('Center').epochList,recType);
            responses(ct,1) = tempCenter.(metric).mean;

            tempSurround = getResponseAmplitudeStats(locationNode.childBySplitValue('Surround').epochList,recType);
            responses(ct,2) = tempSurround.(metric).mean;

            tempCS = getResponseAmplitudeStats(locationNode.childBySplitValue('Center-Surround').epochList,recType);
            responses(ct,3) = tempCS.(metric).mean;
            
            responses(ct,4) = tempCenter.(metric).mean + tempSurround.(metric).mean;
            
            if locationNode.custom.get('isExample')
                addLineToAxis(centerResponse.timeVector,centerResponse.mean,...
                    ['Center_i',num2str(ii),'_l',num2str(ll)],fig1,colors(1,:),'-','none')
                
                addLineToAxis(surroundResponse.timeVector,surroundResponse.mean,...
                    ['Surround_i',num2str(ii),'_l',num2str(ll)],fig2,colors(2,:),'-','none')
                
                addLineToAxis(centerSurroundResponse.timeVector,centerSurroundResponse.mean,...
                    ['CS_i',num2str(ii),'_l',num2str(ll)],fig3,colors(3,:),'-','none')
                
                addLineToAxis(centerResponse.timeVector,linearSumResponse,...
                    ['linSum_i',num2str(ii),'_l',num2str(ll)],fig3,'k','-','none')
                
                addLineToAxis(responses(ct,4),responses(ct,3),...
                    ['egPt',num2str(ii),'_l',num2str(ll)],fig4,colors(3,:),'none','.')
            end
        end
            
        addLineToAxis(responses(:,4), responses(:,3),...
            ['data_i_',num2str(ii)],fig4,'k','none','o')
    end
    limUp = max([responses(:,4)', responses(:,3)']);
    limDown = min([responses(:,4)', responses(:,3)']);
    addLineToAxis([limDown limUp],[limDown limUp],...
        'unity',fig4,'k','--','none')
    
    if ~isempty(figureID)
        makeAxisStruct(fig1,['CSadd_C_', figureID] ,'RFSurroundFigs')

        makeAxisStruct(fig2,['CSadd_S_', figureID] ,'RFSurroundFigs')

        makeAxisStruct(fig3,['CSadd_CS_', figureID] ,'RFSurroundFigs')

        makeAxisStruct(fig4,['CSadd_sum_', figureID] ,'RFSurroundFigs')
    end
    
end
 

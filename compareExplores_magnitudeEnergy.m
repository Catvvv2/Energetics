function fH = compareExplores_magnitudeEnergy(data,behOI,stateOI,norm)

    if nargin < 4
        norm = false;
    end

    global params

%%
    keepM = [];

    fH = figure('Position',[476   608   286   258]); axis square; hold on;

    for f = 1:length(data)
        ores = find(data(f).(stateOI));
        reversals = find(data(f).reversal);

        postWindow = repmat(reversals,1,params.reversalWindow) + repmat([0:params.reversalWindow-1],length(reversals),1);
        postWindow = postWindow(:);

        postOres = intersect(ores,postWindow);
        otherOres = setdiff(ores,postOres);

        tmp = data(f).(behOI);
        if norm; tmp = (tmp - nanmean(tmp))./nanstd(tmp); end

        m = [nanmean(tmp(postOres)),nanmean(tmp(otherOres))];
        e = [nanste(tmp(postOres)),nanste(tmp(otherOres))];
        keepM = [keepM; m];

        h = ttest2(tmp(postOres),tmp(otherOres));

        lH = line([m(1)-e(1),m(1)+e(1)],[m(2),m(2)]);
        lH(2) = line([m(1),m(1)],[m(2)-e(2),m(2)+e(2)]);
        set(lH,'Color',params.colorMap(f,:));
        pH = plot(m(1), m(2),'o',...
            'Color',params.colorMap(f,:),...
            'LineWidth',2,'MarkerSize',10);

        try
        if h
            set(pH,'MarkerFaceColor',params.colorMap(f,:))
        else
            set(pH,'MarkerFaceColor','w')
        end

        catch 
            keyboard
        end
    end

    tmpX = xlim; tmpY = ylim;
    newMin = min([tmpX,tmpY]);
    newMax = max([tmpX,tmpY]);
    ylim([newMin,newMax])
    xlim([newMin,newMax])
    lH = line([newMin newMax],[newMin newMax]);
    set(lH,'Color','k')

    set(gca,'FontSize',16)
    xlabel(['post-change ',stateOI])
    ylabel(['other ',stateOI])

    [p,~,stat] = signrank(diff(keepM'))
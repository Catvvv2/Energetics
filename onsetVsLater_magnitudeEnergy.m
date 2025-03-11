function fH = onsetVsLater_magnitudeEnergy(data,behOI,norm)

    if nargin < 3
        norm = false;
    end

    global params

    %%

    keepM = [];
    fH = figure('Position',[476   608   286   258]); axis square; hold on;

    for f = 1:length(data)
        ores = find(data(f).explore);
        onsets = find(data(f).onsets);
        otherOres = setdiff(ores,onsets);
        nones = setdiff([1:length(data(f).explore)],ores);

        tmp = data(f).(behOI);
        if norm; tmp = (tmp - nanmean(tmp))./nanstd(tmp); end

%         m = [nanmean(tmp(nones)),nanmean(tmp(onsets)),nanmean(tmp(otherOres))];
%         e = [nanmean(tmp(nones)),nanmean(tmp(onsets)),nanmean(tmp(otherOres))];
%         eH = errorbar([1:3],m,e);
%         set(eH,'Color',params.colorMap(f,:))

        h = ttest2(tmp(onsets),tmp(otherOres));

        m = [nanmean(tmp(onsets)),nanmean(tmp(otherOres))];
        e = [nanste(tmp(onsets)),nanste(tmp(otherOres))];
        keepM = [keepM; m];

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
    xlabel(['onsets'])
    ylabel(['other explores'])


    [p,~,stat] = signrank(diff(keepM'))
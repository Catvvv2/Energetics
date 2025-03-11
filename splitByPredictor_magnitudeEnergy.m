function fH = splitByPredictor_magnitudeEnergy(data,behOI,predOI,z_score)
    % looks as whether there's a difference in behOI, depending on the
    % predictor variable predOI (must be binary)

    global params

    fH = figure(); hold all;
    
    M = []; sigSessions = 0;

    for f = 1:length(data)
    
        y = data(f).(behOI);
    
        if nargin > 3 && z_score
            y = (y - nanmean(y))./nanstd(y);
        end
    
        m = [nanmean(y(data(f).(predOI)==0)),...
            nanmean(y(data(f).(predOI)==1))];
        e = [nanste(y(data(f).(predOI)==0)),...
            nanste(y(data(f).(predOI)==1))];
    
        if params.plotRaw
        
            h = errorbar([0,1],m,e);
            set(h,'Color',params.colorMap(f,:));
            if ttest2(y(data(f).(predOI)==0),y(data(f).(predOI)==1))
                plot(1.2,m(2),'*','Color',params.colorMap(f,:))
                sigSessions = sigSessions + 1;
            end
        end

        M = [M; m];
    
    end
    
    % spit some information about the analysis to the command window
    sprintf('n. sig. sessions = %2.0f / %2.0f',sigSessions,length(data))

    % append axis labels
    set(gca,'FontSize',params.FontSize)
    xlabel(predOI);
    if nargin > 3 && z_score
        ylabel(strcat(behOI,' (z-scored)'))
    else
        ylabel(behOI);
    end
    set(gca,'XTick',[0,1]);
    xlim([-0.5 1.5]);

    h = errorbar([0,1],nanmean(M),nanste(M));
    set(h,'LineWidth',2,'Color',params.defaultBlack)

    effectSize = nanmean(M(:,1)-M(:,2))
    effectRange = [min(M(:,1)-M(:,2)),max(M(:,1)-M(:,2))]
    % [h,p,ci,stats] = ttest(M(:,1),M(:,2))
    [p,h,stats] = signrank(M(:,1),M(:,2))
    plot(0.5,nanmean(nanmean(M))+0.1,'*k')
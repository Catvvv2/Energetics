function fH = blockChange_magnitudeEnergy(data,behOI,z_score,blockChangeType)
    % takes behOI string and plots that field in data struct aligned to the
    % block changes
    
    if nargin < 4
        blockChangeType = 'blockChange';
    end

    if nargin < 3
        z_score = false;
    end

    global params

    fH = figure(); hold all;
    
    % export parameters
    minTrial = params.minTrial;
    maxTrial = params.maxTrial;

    xpos = minTrial:maxTrial;
    keepY = []; sessIdx = []; ref = [];
    
    for f = 1:length(data)
    
        % pull out the y data
        y = data(f).(behOI);
        if z_score
            y = (y - nanmean(y))./nanstd(y);
        end

        % aligned to objective block changes
        changes = find([data(f).(blockChangeType)]==1); changes = changes(:);
    
        % exclude too longs/too shorts
        changes = changes(and(changes+minTrial > 0,changes+maxTrial < length([data(f).reversal])));    
    
        % indexes now
        indx = repmat(changes,1,size(xpos,2)) + repmat(xpos,size(changes,1),1);
    
        if params.plotRaw
            % now the plot of this session's average
            plot(xpos,nanmean(y(indx)),'Color',params.colorMap(f,:));
        end

        keepY = [keepY; y(indx)];
        ref = [ref; nanmean(y)];
        sessIdx = [sessIdx; repmat(f,size(y(indx),1),1)];
    end

    % spit some information about the analysis to the command window
    sprintf('errorbars = across trial sequences \n N = %2.0f',size(keepY,1))

    % now we'll do our ribbon plots
    m = nanmean(keepY);
    e = nanste(keepY);

    plot(xpos,m,'Color',params.defaultBlack,...
        'LineWidth',2)
    plot(xpos,m+e,'--','Color',params.defaultBlack,...
        'LineWidth',1)
    plot(xpos,m-e,'--','Color',params.defaultBlack,...
        'LineWidth',1)
 
    % append axis labels
    set(gca,'FontSize',params.FontSize)
    if nargin > 2 && z_score
        ylabel(strcat(behOI,' (z-scored)'))
        line([xpos(1), xpos(end)],[0 0],'Color','k')
    else
        ylabel(behOI);
    end
    xlabel(['trials since ',blockChangeType])

    % reference line at reversal trial
    y_tmp = ylim;
    line([0 0],[y_tmp(1) y_tmp(2)],'Color','k')

    % now we want to do some significance testing
    % we'll compare each post-change point w/ the baseline avg
    reference = mean(keepY(:,xpos < 0),2);
%     reference = mean(keepY(:,xpos < 0),2);
    reference = mean(ref);
%     keyboard()

    for t = find(xpos>=-5)
        h = ttest(keepY(:,t),reference);
        if h
            plot(xpos(t),y_tmp(2),'*k')
        end
    end

    
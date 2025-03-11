function [fH,data] = spikingPlots_magnitudeEnergy(data,predOI)

    global params

    fH = figure();
    
    M = [];
    sparsity = NaN(length(data),2);
    cmap = get(gca,'ColorOrder'); cmap = [cmap;cmap];
    
    for f = 1:length(data)
    
        spkRate = data(f).spikeRate;
        spkRate = (spkRate-repmat(min(spkRate,[],2),1,size(spkRate,2))) ./ repmat(max(spkRate,[],2)-min(spkRate,[],2),1,size(spkRate,2));
    
        % now for some plots, first the single neurons
        subplot(2,3,[1,2,4,5]); axis square; hold all;
        colororder(params.colorMap);
    
        ph = plot(nanmean(spkRate(:,data(f).(predOI)==0),2),...
            nanmean(spkRate(:,data(f).(predOI)==1),2),'.');
    
        M = [M; nanmean(spkRate(:,data(f).(predOI)==0),2),...
            nanmean(spkRate(:,data(f).(predOI)==1),2)];
    
        h = ttest(nanmean(spkRate(:,data(f).(predOI)==1),2),nanmean(spkRate(:,data(f).(predOI)==0),2));
        if h
            plot(0.2+f/50,0.9,'*','Color',get(ph,'Color'))
        end

        % now sparsity
        sparsity(f,1) = nanmean(data(f).trialSparsity(data(f).(predOI)==0));
        sparsity(f,2) = nanmean(data(f).trialSparsity(data(f).(predOI)==1));
    
        m = sparsity(f,:);
        e = [nanste(data(f).trialSparsity(data(f).(predOI)==0)),...
            nanste(data(f).trialSparsity(data(f).(predOI)==0))];
    
        subplot(2,3,6); axis square; hold all;
        colororder(params.colorMap);
        h = errorbar([0,1],m,e,'Color',cmap(f,:));
        
        if ttest2(data(f).trialSparsity(data(f).(predOI)==0),data(f).trialSparsity(data(f).(predOI)==1))
            plot(1.2,m(2),'*','Color',get(h,'Color'))
        end
    
    end
    
    subplot(2,3,[1,2,4,5]);
    set(gca,'FontSize',params.FontSize)
    line([0, 1],[0, 1]);
    xlabel([predOI,'== 0']); ylabel([predOI,'== 1']);

    subplot(2,3,3); axis square; hold on;
    set(gca,'FontSize',params.FontSize)
    tmp = M(:,1) - M(:,2);
    h = histogram(tmp,'normalization','probability');
    set(h,'LineStyle','none','FaceColor',[.5 .5 .5])
    h = line([0 0],[0 0.15]); set(h,'Color','k');
    plot(nanmean(tmp))
    xlabel({'difference in spike rate';[predOI,' == 0 - 1']})
    [h,p] = ttest(tmp);
    title(sprintf('p < %0.4f',p))
    ylabel('p(neuron)')
    xlim([-0.15 0.15])
    
    subplot(2,3,6);
    set(gca,'FontSize',params.FontSize)
    xlim([-0.5 1.5]);
    ylabel('sparsity')
    xlabel(predOI);
    set(gca,'XTick',[0,1]);

    effectSize = mean(diff(sparsity'))
    effectRange = [min(diff(sparsity')),max(diff(sparsity'))]
    [p,h,stats] = signrank(diff(sparsity'))
    title(sprintf('p < %0.2f',ceil(p * 100)/100))


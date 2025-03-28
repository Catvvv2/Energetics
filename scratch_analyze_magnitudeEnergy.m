%% load the datafile

cd('/Users/becket/Documents/MATLAB/magnitudeEnergy')
clear;
analysisParams_magnitudeEnergy;

cd(params.dataDir);
load(params.dataFile);


%% messing w/ the session-by-session measures, re-cacluating

for f = 1:length(data)
%     % quickly re-calculate to be exclusive
%     data(f).sparsityDifference = nanmean(data(f).trialSparsity(data(f).explore==1),2)-nanmean(data(f).trialSparsity(data(f).explore==0),2);
% 
%     m = [nanmean(data(f).trialSpikeRate(data(f).explore==0),2),...
%         nanmean(data(f).trialSpikeRate(and(data(f).explore==1,data(f).onsets==0)),2)]
%     data(f).spikeRateDifference = diff(m);
%     
%     m = [nanmean(data(f).magnitude(and(data(f).explore==0,data(f).onsets==0))),...
%         nanmean(data(f).magnitude(and(data(f).explore==1,data(f).onsets==0)))];
% 
%     m = [nanmean(data(f).magnitude(data(f).explore==0)),...
%         nanmean(data(f).magnitude(and(data(f).explore==1,data(f).onsets==0)))];
%     data(f).magnitudeDifference = diff(m) ./ std(data(f).magnitude); % ore vs oit

    % add offset labels
    data(f).offsets = [0 diff(data(f).stateLabels == 1)] == -1;
    % if you look at the timeseries - the offset peaks on the trial before

    m = [nanmean(data(f).magnitude(data(f).offsets==0)),...
        nanmean(data(f).magnitude(data(f).offsets==1))];
    data(f).magnitudeDifference_offset = diff(m) ./ std(data(f).magnitude); % ore vs oit

    m = [nanmean(data(f).trialSpikeRate(data(f).offsets==0)),...
        nanmean(data(f).trialSpikeRate(data(f).offsets==1))];
    data(f).spikeRateDifference_offset = diff(m) ./ std(data(f).magnitude); % ore vs oit


end

%%
yStr = 'magnitude'; ylimits = [-0.1 1.25];
% yStr = 'spikeRate'; ylimits = [0 0.025];

% yStr = 'sparsity';

figure('Position',[476   453   891   413]);
for plt = 1:6
    switch plt
        case 1
            predVar = 'basinDifference';
            xlimits = [1.1 2];
            dVar = [yStr,'Difference'];
            tStr = 'explore';
        case 2
            predVar = 'basinDifference';
            xlimits = [1.1 2];
            dVar = [yStr,'Difference_offset'];
            tStr = 'offset'; 
        case 4
            predVar = 'barrierHeight';
            xlimits = [1.3 2.4];
            dVar = [yStr,'Difference'];
            tStr = 'explore'; 
        case 5
            predVar = 'barrierHeight';
            xlimits = [1.3 2.4];
            dVar = [yStr,'Difference_offset'];
            tStr = 'offset';
        case 3
            predVar = 'basinDifference';
            xlimits = [1.1 2];
            dVar = [yStr,'Difference_reward'];
            tStr = 'reward';
        case 6
            predVar = 'barrierHeight';
            xlimits = [1.3 2.4];
            dVar = [yStr,'Difference_reward'];
            tStr = 'reward';
    end

    if strcmp(predVar,'barrierHeight'); xlimits = [1.3 2.4];
    elseif strcmp(predVar,'basinDifference'); xlimits = [1.1 2];
    else; xlimits = [1 3]; end

    subplot(2,3,plt); hold all;
    xlim(xlimits);
%     tmp = ylim; ylim([tmp(1)-0.1 tmp(end)+0.1]);
    
    y = [data.(dVar)];
    if strcmp(tStr,'reward'); y = 1-y; end % recode to omissions
    ph = plot([data.(predVar)],y,...
        '.k','MarkerSize',20)
    lsline;

    % or we can color the dots
    delete(ph);
    for f = 1:length(data)
        plot([data(f).(predVar)],y(f),...
            '.','MarkerSize',20,'Color',params.colorMap(f,:))
    end
    
    xlabel(predVar)
    ylabel({['difference in ',yStr],'1 minus 0'})
    
    [r,p] = corr([data.(predVar)]',[data.(dVar)]','Type','Spearman');
    title([tStr, sprintf(', rho = %0.2f, p < %0.4f',r,p)]);

    
end

[r,p] = corr([data.basinDifference]',[data.barrierHeight]')

[r,p] = corrcoef([data.magnitudeDifference_onset],[data.magnitudeDifference_reward])

print([params.figDir,'/',yStr,'CorrWithinSession','_analyze_magnitudeEnergy'],params.figFormat)




%% now looking at the 1001 person dataset, to see if the barrier height is related to the reward history filter

% fH = barrierHeightByRwdHistoryFilter_humanManifolds;

% as long as we're here, we can quickly calculate the effect of reward on
% the p(onset) for this session; we'll do this as a simple index:
for k = 1:length(data)
    data(k).rwdFilterA = (nanmean(data(k).onsets(data(k).last_rwd==0))-nanmean(data(k).onsets(data(k).last_rwd==1))) ./ ...
        nanmean(data(k).onsets);
end
% session by session effects of rewards on onsets


figure(); 
% append the data from the monkeys here
hold on;
for analysis = 1:4
    switch analysis
        case 1
            xStr = 'barrierHeight';
            y = [data.rwdFilterA];
            tStr = 'rwd sensitivity';
        case 2
            xStr = 'basinDifference';
            y = [data.rwdFilterA];
        case 3
            xStr = 'barrierHeight';
            y = [data.rwdFilterB];
            tStr = 'exp param';
        case 4
            xStr = 'basinDifference';
            y = [data.rwdFilterB];
    end

    x = [data.(xStr)];

    subplot(2,2,analysis); hold on;
    pH = plot(x,y,'.k'); lsline; delete(pH)

    for f = 1:length(data)
        plot(x(f),y(f),...
            '.','MarkerSize',20,'Color',params.colorMap(f,:))
    end

    [r,p] = corr(x',y','Type','Spearman');
    title([tStr, sprintf(', rho = %0.2f, p < %0.4f',r,p)]);

    xlabel(xStr);

    % THIS IS NOT WORKING FOR SOME REASON????
    if strcmp(xStr,'barrierHeight'); xlimits = [1.3 2.4];
    elseif strcmp(xStr,'basinDifference'); xlimits = [1.1 2];
    else; xlimits = [1 3]; end

end

% % there is a systematic relationship between the reward history filter and
% % the barrier height/basin difference - it's super clear in the 1001 person
% % dataset, though there's some issues w/ fitting there too -- need a better
% % way to verify that the model is actually fitting ok
% 
% % [r,p] = corrcoef([data.kernelB],[data.spikeRateDifference_onset])
% % [r,p] = corrcoef([data.kernelB],[data.magnitudeDifference_onset])
% % 
% % [r,p] = corrcoef([data.kernelLast],[data.spikeRateDifference_onset])
% % [r,p] = corrcoef([data.kernelLast],[data.magnitudeDifference_onset])
% 
% % [r,p] = corr([data.kernelA]',[data.barrierHeight]','Type','Spearman')
% % [r,p] = corr([data.kernelLast]',[data.barrierHeight]','Type','Spearman')
% [r,p] = corr([data.rwdFilterB]',[data.barrierHeight]','Type','Spearman')
% 
% % [r,p] = corr([data.kernelA]',[data.basinDifference]','Type','Spearman')
% % [r,p] = corr([data.kernelLast]',[data.basinDifference]','Type','Spearman')
% [r,p] = corr([data.rwdFilterB]',[data.basinDifference]','Type','Spearman')

%%

fNames = {'barrierHeight',...
    'magnitudeDifference_reward',...
    'rwdFilterA'};

X = [];
for k = 1:length(fNames);
    X = [X, [data.(fNames{k})]'];
end

[r,p] = partialcorr(X,'type','Pearson')

% [r,p] = corr(X,'type','Spearman')


%% append offset labels?
% looks like any appearance of an offset transient is just due to
% exploration happening on the previous trial!!

for k = 1:length(data)
    data(k).offsets = [0 diff(data(k).stateLabels == 1)] == -1;
end

for analysis = 1:4
    switch analysis
        case 1
            behOI = 'explore'; norm = true;
        case 2
            behOI = 'explore'; norm = true;
        case 3
            behOI = 'trialSpikeRate'; norm = true;
        case 4
            behOI = 'trialSparsity'; norm = true;
    end
    fH = blockChange_magnitudeEnergy(data,behOI,norm,'offsets');
    
    set(fH,'Position',[476   632   371   234])
    print([params.figDir,'/',behOI,'OffsetAlignedByTrial','_analyze_magnitudeEnergy'],params.figFormat)

end
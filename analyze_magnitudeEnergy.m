%% load the datafile

cd('/Users/becket/Documents/MATLAB/magnitudeEnergy')
clear;
analysisParams_magnitudeEnergy;

cd(params.dataDir);
load(params.dataFile);

%% append variables

% trial spike variance is avg distance from avg. firing rate w/in each neuron

for k = 1:length(data)
    data(k).trialSpikeVariance = nanmean((data(k).spikeRate - repmat(nanmean(data(k).spikeRate,2),1,size(data(k).spikeRate,2))).^2);
end

%% illustrate the HMM model fits

fH = baselinePlots_magnitudeEnergy(data);
print([params.figDir,'/','baselineHMMresults_analyze_magnitudeEnergy'],params.figFormat)

%% reward maximizing vs information maximizing?

[fH1,fH2] = rewardVsInfoMax_magnitudeEnergy(data);
figure(fH1);
print([params.figDir,'/','rewardVsInfoMax_analyze_magnitudeEnergy'],params.figFormat)
figure(fH2);
print([params.figDir,'/','rewardVsInfoMaxChUnCh_analyze_magnitudeEnergy'],params.figFormat)


%% block-change aligned analyses, behavior

for analysis =  1:4
    switch analysis
        case 1
            behOI = 'reward'; norm = false;
        case 2
            behOI = 'choseBest'; norm = false;
        case 3
            behOI = 'explore'; norm = false;
        case 4
            behOI = 'onsets'; norm = false;
    end

    fH = blockChange_magnitudeEnergy(data,behOI,norm);
    set(fH,'Position',[476   632   371   234])

    print([params.figDir,'/',behOI,'ByBlock_analyze_magnitudeEnergy'],params.figFormat)
end

%% curious about differences between the 2 kinds of block changes


tStr = {'reversal','newImages'};

behOI = 'explore'; norm = false; ybits = NaN;

% behOI = 'choseBest'; norm = false; ybits = [0 1];

for analysis = 1:2
    fH = blockChange_magnitudeEnergy(data,behOI,norm,tStr{analysis});
    title(tStr{analysis});
%     ylim([0 0.8])
    set(fH,'Position',[476   632   371   234])

    if ~isnan(ybits)
        ylim([ybits])
    end

    print([params.figDir,'/',behOI,'ByBlock',tStr{analysis},'_analyze_magnitudeEnergy'],params.figFormat)

end

%% we want to understand the relationship between exploration and last reward?

% first, they tend to co-occur
behOI = 'explore';
% predOI = 'last_rwd';
predOI = 'reward';
% predOI = 'choseBest';

fH = splitByPredictor_magnitudeEnergy(data,behOI,predOI);
set(fH,'Position',[476   642   219   224])

print([params.figDir,'/','exploreBy',predOI,'_analyze_magnitudeEnergy'],params.figFormat)


%% single unit analyses
% last reward and exploration both (1) increase

% predOI = 'last_rwd'
predOI = 'explore';

[fH,data] = spikingPlots_magnitudeEnergy(data,predOI) % also appends sparsity, spike rate
print([params.figDir,'/','spikingAnalyses_analyze_magnitudeEnergy'],params.figFormat)

% then statistics to check to see if exploration and reward have additive effects on each
behOI = 'trialSpikeRate';
dissociateStateFromReward_magnitudeEnergy(data,behOI)

% same for sparsity
behOI = 'trialSparsity';
dissociateStateFromReward_magnitudeEnergy(data,behOI)

%% do a session-by-session plot of the same

behOI = 'trialSpikeRate'; 
behOI = 'trialSparsity'; 
behOI = 'trialSpikeVariance';

predOI = 'explore';
fH = splitByPredictor_magnitudeEnergy(data,behOI,predOI,true);
set(fH,'Position',[476   642   219   224])

print([params.figDir,'/',behOI,'ByExplore','_analyze_magnitudeEnergy'],params.figFormat)

%% block aligned plots for neural activity

for analysis =  1:4
    switch analysis
        case 1
            behOI = 'trialSpikeRate'; norm = true;
        case 2
            behOI = 'trialSparsity'; norm = true;
        case 3
            behOI = 'magnitude'; norm = true;
        case 4
            behOI = 'trialSpikeVariance'; norm = true;
    end

    fH = blockChange_magnitudeEnergy(data,behOI,norm);
    set(fH,'Position',[476   632   371   234])

    print([params.figDir,'/',behOI,'ByBlock_analyze_magnitudeEnergy'],params.figFormat)
end

%% sparsity only peaks w/ new images:

behOI = 'trialSparsity'; norm = true;
% behOI = 'trialSpikeRate'; norm = true;
% behOI = 'magnitude'; norm = true;

fH = blockChange_magnitudeEnergy(data,behOI,norm,'reversal');
set(fH,'Position',[476   632   371   234])
print([params.figDir,'/',behOI,'TimeCourse_reversal_analyze_magnitudeEnergy'],params.figFormat)
fH = blockChange_magnitudeEnergy(data,behOI,norm,'newImages');
set(fH,'Position',[476   632   371   234])
print([params.figDir,'/',behOI,'TimeCourse_newImages_analyze_magnitudeEnergy'],params.figFormat)

%% population magnitude analyses

behOI = 'magnitude'; % magnitude itself has additive effects
% behOI = 'magnitude_change'; % change is *only* related to reward

predOI = 'last_rwd';
fH = splitByPredictor_magnitudeEnergy(data,behOI,predOI,true);
set(fH,'Position',[476   642   219   224])

print([params.figDir,'/',behOI,'ByLastRwd','_analyze_magnitudeEnergy'],params.figFormat)

predOI = 'explore';
fH = splitByPredictor_magnitudeEnergy(data,behOI,predOI,true);
set(fH,'Position',[476   642   219   224])

print([params.figDir,'/',behOI,'ByExplore','_analyze_magnitudeEnergy'],params.figFormat)

% estimate relative contribution of reward and state
dissociateStateFromReward_magnitudeEnergy(data,behOI)

%% relationship between spike rate and magnitude

binWidth = .025;
nSamples = 20;

figure('Position',[476   739   181   127]); hold on;

[rs,ps] = deal(NaN(length(data),nSamples));
for f = 1:length(data)
    for s = 1:nSamples
        idx = randi(length(data(f).magnitude),[length(data(f).magnitude),1]);
        [r,p] = corrcoef(data(f).trialSpikeRate(idx),data(f).magnitude(idx));
        rs(f,s) = r(1,2); ps(f,s) = p(1,2);
    end

    hH = histogram(rs(f,:),'BinEdges',[0:binWidth:1]);
    set(hH,'LineStyle','none','FaceColor',params.colorMap(f,:))
end

min(nanmean(rs')')
max(nanmean(rs')')

set(gca,'FontSize',16)
ylabel('counts')
xlabel('correlation')

print([params.figDir,'/','corrBtwMagnitudeAndSpiking','_analyze_magnitudeEnergy'],params.figFormat)


%% determine if there are any systematic differences between post-block change behaviors
%   and behaviors that occurred elsewhere in the block

% stateOI = 'onsets';
% params.reversalWindow = 5; % 

stateOI = 'explore';
params.reversalWindow = 5; %

norm = true;

for analysis = 1:4
    switch analysis
        case 1
            behOI = 'last_rwd';
        case 2
            behOI = 'trialSparsity';
        case 3
            behOI = 'trialSpikeRate';
        case 4
            behOI = 'magnitude';
    end

    behOI

    fH = compareExplores_magnitudeEnergy(data,behOI,stateOI,norm);
    title(behOI);

    print([params.figDir,'/',behOI,'Compare',stateOI,'_analyze_magnitudeEnergy'],params.figFormat)
end

%% now focus specifically on the onsets

predOI = 'onsets';
fH = splitByPredictor_magnitudeEnergy(data,behOI,predOI,true);
set(fH,'Position',[476   642   219   224])

print([params.figDir,'/',behOI,'ByOnset','_analyze_magnitudeEnergy'],params.figFormat)

% estimate relative contribution of reward and state
dissociateStateFromReward_magnitudeEnergy(data,behOI,predOI)

%% are the onsets more high FR than the later trials?

norm = true;

for analysis = 1:3
    switch analysis
        case 1
            behOI = 'trialSparsity';
        case 2
            behOI = 'trialSpikeRate';
        case 3
            behOI = 'magnitude';
%             behOI = 'last_rwd';
%             behOI = 'trialSpikeVariance';
    end
    behOI

    fH = onsetVsLater_magnitudeEnergy(data,behOI,norm);
    title(behOI);

    print([params.figDir,'/',behOI,'onsetVsLater_analyze_magnitudeEnergy'],params.figFormat)
end

%% now control for reward effects
% they are less likely to be rewarded before onsets than other explores
%    so reward effects could just be driving the onset difference

behOI = 'magnitude';
dissociateOnsetExploreAndReward_magnitudeEnergy(data,behOI)

behOI = 'trialSpikeRate';
dissociateOnsetExploreAndReward_magnitudeEnergy(data,behOI)
% magnitude, doesn't survive, but trial spike rate does

%% aligned to the onsets?

for analysis = 1:3
    switch analysis
        case 1
            behOI = 'magnitude'; norm = true;
        case 2
            behOI = 'trialSpikeRate'; norm = true;
        case 3
            behOI = 'trialSparsity'; norm = true;
    end
    fH = blockChange_magnitudeEnergy(data,behOI,norm,'onsets');
    
    set(fH,'Position',[476   632   371   234])
    print([params.figDir,'/',behOI,'OnsetAlignedByTrial','_analyze_magnitudeEnergy'],params.figFormat)

end

%% bump graphs for each session??

params.bumpA = 1; % scaling constant
params.bumpWrtOre = false; % else oit; % must be true, w.r.t. ORE, for barrier height to be sig

figure('Position',[476   573   560   293]);
subplot(2,3,[1:2,4:5]);hold all;

Eout = NaN(length(data),3);
for k = 1:length(data)
    T = data(k).T;
    pi = data(k).stationaryDist;
    if params.bumpWrtOre
        [E,~,E_b] = twoWellPotential_v2([pi(1),1-pi(1)],[1-T(1,1),T(2,1)],false,params.bumpA); % w.r.t. ORE
    else
        [E,~,E_b] = twoWellPotential_v2([1-pi(1),pi(1)],[1-T(2,2),1-T(1,1)],false,params.bumpA); % w.r.t to OIT
    end

    % 3nd part of this normalizes to the mean
    Eout(k,:) = E([1,2,end]) - nanmean(E([1,2,end]));
    
    plot(Eout(k,:),'Color',params.colorMap(k,:))
end

xlim([0.5 3.5]);

subplot(2,3,[3]);hold all;
for k = 1:length(data)
    bar(k,data(k).basinDifference,...
        'LineStyle','none','FaceColor',params.colorMap(k,:));
end
ylim([-3,3])
set(gca,'XTick',[1:length(data)],'FontSize',12)
xlabel('session #'); ylabel('\Delta energy')

subplot(2,3,[6]);hold all;
for k = 1:length(data)
    bar(k,data(k).barrierHeight,...
        'LineStyle','none','FaceColor',params.colorMap(k,:));
end
ylim([-3,3])
set(gca,'XTick',[1:length(data)],'FontSize',12)
xlabel('session #'); ylabel('activation energy')

print([params.figDir,'/','bumpGraphs','_analyze_magnitudeEnergy'],params.figFormat)


%% session-by-session relationship between thermodynamic measures and the neural activity

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
            dVar = [yStr,'Difference_onset'];
            tStr = 'onset'; 
        case 4
            predVar = 'barrierHeight';
            xlimits = [1.3 2.4];
            dVar = [yStr,'Difference'];
            tStr = 'explore'; 
        case 5
            predVar = 'barrierHeight';
            xlimits = [1.3 2.4];
            dVar = [yStr,'Difference_onset'];
            tStr = 'onset';
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


%% do rewards do work on the system? we'll try to get at this w/ a mediation analysis
% also, Wiener filter; a temporal relationship: omitted rewards proceed onsets

params.nLagsWiener = 10;
params.nParamsWiener = 1;

[fH,data] = rewardsWork_magnitudeEnergy(data);

print([params.figDir,'/','rwdHistoryFilter','_analyze_magnitudeEnergy'],params.figFormat)

function [fH1,fH2] = rewardVsInfoMax_magnitudeEnergy(data);

global params

% ok, so the idea is that we want to know what the information gain is for
% each choice, given some previous sequence of choices
% we'll do this very simply, but just making a Bayesian reward-tracking
% model for each feature dimension

nOpts = 2; evOut = []; hOut = []; infoOut = [];
[fH1] = figure(99); set(gcf,'Position',[476   574   742   292]);

for k = 1:length(data)

    newImages = [data(k).newImages];
    imChoice = [data(k).imChoice]+1;
    locChoice = [data(k).locChoice]+1;
    rewards = [data(k).reward];
    
    [im_alpha,im_beta,im_EV,...
        loc_alpha,loc_beta,loc_EV] = deal(NaN(length(newImages)+1,nOpts));
    
    for t = 1:length(newImages)
        if newImages(t) % restart prior when the images switch
            im_alpha(t,:) = 1; im_beta(t,:) = 1;
            loc_alpha(t,:) = 1; loc_beta(t,:) = 1;
        end
    
        % carry forward unchosen values
        im_alpha(t+1,:) = im_alpha(t,:); im_beta(t+1,:) = im_beta(t,:);
        loc_alpha(t+1,:) = loc_alpha(t,:); loc_beta(t+1,:) = loc_beta(t,:);
    
        % when chosen, increment alpha or beta, depending on reward
        im_alpha(t+1,imChoice(t)) = im_alpha(t,imChoice(t)) + double(rewards(t) > 0);
        im_beta(t+1,imChoice(t)) = im_beta(t,imChoice(t)) + double(rewards(t) == 0);
        loc_alpha(t+1,imChoice(t)) = loc_alpha(t,imChoice(t)) + double(rewards(t) > 0);
        loc_beta(t+1,imChoice(t)) = loc_beta(t,imChoice(t)) + double(rewards(t) == 0);
    
    end
    
    % calculate the expected value
    im_EV = im_alpha ./ (im_alpha + im_beta);
    loc_EV = loc_alpha ./ (loc_alpha + loc_beta);
    EV = im_EV + loc_EV;
    
    % we also want to calculate info gain
    % as the E(reduction in entropy) from making a choice
    [chosenEV,chosenH,chosenInfo,...
        unchosenEV,unchosenH,unchosenInfo] = deal(zeros(length(im_alpha)-1,1));
    for dim = 1:2
        switch dim
            case 1
                alpha = im_alpha; beta = im_beta; choice = imChoice;
            case 2
                alpha = loc_alpha; beta = loc_beta; choice = locChoice;
        end
    
        EV = alpha ./ (alpha + beta);
        H = -(EV.*log2(EV) + (1-EV).*log2((1-EV)));
    
        posEV = (alpha+1) ./ ((alpha+1) + beta);
        posH = - (posEV.*log2(posEV) + (1-posEV).*log2(1-posEV));
        negEV = alpha ./ (alpha + (beta+1));
        negH = - (negEV.*log2(negEV) + (1-negEV).*log2(1-negEV));
        futureH = EV.*posH + (1-EV).*negH;
        Info = H - futureH;

        % we'll assume additivity in the 2 dimensions
        chosenIndex = sub2ind(size(EV),1:length(choice),choice');
        unchosenIndex = sub2ind(size(EV),1:length(choice),mod(choice',2)+1);
        chosenEV = chosenEV+EV(chosenIndex)';
        chosenH = chosenH+H(chosenIndex)';
        chosenInfo = chosenInfo+Info(chosenIndex)';

        unchosenEV = unchosenEV+EV(unchosenIndex)';
        unchosenH = unchosenH+H(unchosenIndex)';
        unchosenInfo = unchosenInfo+Info(unchosenIndex)';

        % we could calculate an index of chosen / unchosen here, but will
        % do this later instead
%         chosenEV = chosenEV+[(EV(chosenIndex)-EV(unchosenIndex))./mean(EV,2)]';
%         chosenH = chosenH+[(H(chosenIndex)-H(unchosenIndex))./mean(H,2)]';
%         chosenInfo = chosenInfo+[(Info(chosenIndex)-Info(unchosenIndex))./mean(Info,2)]';
    end

    % only for EV, we will calculate the expected value
    chosenEV = chosenEV ./ 2;
    unchosenEV = unchosenEV ./ 2;
    % info theoretic measures, we'll assume are additive

    % save stuff to the datafile while we're here
    data(k).chosenEV = chosenEV;
    data(k).chosenH = chosenH;
    data(k).chosenInfo = chosenInfo;
    data(k).unchosenEV = unchosenEV;
    data(k).unchosenH = unchosenH;
    data(k).unchosenInfo = unchosenInfo;
    
    % now plots
    for plt = 1:3
        figure(99);
        switch plt
            case 1
                behOI = chosenEV;
            case 2
                behOI = chosenH;
            case 3
                behOI = chosenInfo;
        end
    
        subplot(1,3,plt); hold on; axis square;
        m = [nanmean(behOI(data(k).explore==0)),...
            nanmean(behOI(data(k).explore==1))];
        e = [nanste(behOI(data(k).explore==0)),...
            nanste(behOI(data(k).explore==1))];

        markerStr = '.'; markerSize = 25;
        plot(m(1),m(2),markerStr,'Color',params.colorMap(k,:),...
            'MarkerSize',markerSize)
        plot([m(1),m(1)],[m(2)-e(2),m(2)+e(2)],'-k')
        plot([m(1)-e(1),m(1)+e(1)],[m(2),m(2)],'-k')

        switch plt
            case 1
                evOut = [evOut;m];
            case 2
                hOut = [hOut;m];
            case 3
                infoOut = [infoOut;m];
        end
    end

end

% clean up the plots
figure(99);
tStr = {'EV','H','Info'};
for plt = 1:3
    subplot(1,3,plt);
    tmpY = ylim; tmpX = xlim;
    xlim([min([tmpY(1),tmpX(1)]),max([tmpY(2),tmpX(2)])])
    ylim([min([tmpY(1),tmpX(1)]),max([tmpY(2),tmpX(2)])])
    line([min([tmpY(1),tmpX(1)]),max([tmpY(2),tmpX(2)])],...
        [min([tmpY(1),tmpX(1)]),max([tmpY(2),tmpX(2)])],'Color','k')
    set(gca,'FontSize',12)
    title(tStr{plt});
    xlabel('rule state')
    ylabel('explore state')
end

% spit out the stats, this is to tell if there's more of this variable
% during exploit than explore
for analysis = 1:3
    switch analysis
        case 1
            out = evOut; behStr = 'EV';
        case 2
            out = hOut; behStr = 'H';
        case 3
            out = infoOut; behStr = 'Info';
    end

    [p,h,stat] = signrank(out(:,1)-out(:,2));
    sprintf([behStr,' oit-ore diff p = %0.4f, sign = %2.0f'],p,stat.signedrank)
    sprintf('mean diff = %2.4f; range = [%2.4f, %2.4f]',mean(out(:,1)-out(:,2)),...
        min(out(:,1)-out(:,2)),max(out(:,1)-out(:,2)))

end

%% now compare chosen to unchosen

[fH2] = figure(98); set(gcf,'Position',[476   574   742   272]);
[oitM,oreM] = deal(NaN(length(data),2,3));

for k = 1:length(data)

    for plt = 1:3
        switch plt
            case 1
                chosen = data(k).chosenEV;
                unchosen = data(k).unchosenEV;
            case 2
                chosen = data(k).chosenH;
                unchosen = data(k).unchosenH;
            case 3
                chosen = data(k).chosenInfo;
                unchosen = data(k).unchosenInfo;
        end
    
        subplot(2,3,plt); hold on; axis square;
        m = [nanmean(chosen(data(k).explore==0)),...
            nanmean(unchosen(data(k).explore==0))];
        e = [nanste(chosen(data(k).explore==0)),...
            nanste(unchosen(data(k).explore==0))];
    
        markerStr = '.'; markerSize = 25;
        plot(m(1),m(2),markerStr,'Color',params.colorMap(k,:),...
            'MarkerSize',markerSize)
        plot([m(1),m(1)],[m(2)-e(2),m(2)+e(2)],'-k')
        plot([m(1)-e(1),m(1)+e(1)],[m(2),m(2)],'-k')

        oitM(k,:,plt) = m;

        subplot(2,3,plt+3); hold on; axis square;
        m = [nanmean(chosen(data(k).explore==1)),...
            nanmean(unchosen(data(k).explore==1))];
        e = [nanste(chosen(data(k).explore==1)),...
            nanste(unchosen(data(k).explore==1))];
    
        plot(m(1),m(2),markerStr,'Color',params.colorMap(k,:),...
            'MarkerSize',markerSize)
        plot([m(1),m(1)],[m(2)-e(2),m(2)+e(2)],'-k')
        plot([m(1)-e(1),m(1)+e(1)],[m(2),m(2)],'-k')

        oreM(k,:,plt) = m;
    end
end

for plt = 1:6
    subplot(2,3,plt);
    tmpY = ylim; tmpX = xlim;
    xlim([min([tmpY(1),tmpX(1)]),max([tmpY(2),tmpX(2)])])
    ylim([min([tmpY(1),tmpX(1)]),max([tmpY(2),tmpX(2)])])
    line([min([tmpY(1),tmpX(1)]),max([tmpY(2),tmpX(2)])],...
        [min([tmpY(1),tmpX(1)]),max([tmpY(2),tmpX(2)])],'Color','k')
    set(gca,'FontSize',12)
    if plt > 3; xlabel('chosen'); end
    ylabel('unchosen')
end

for analysis = 1:3

    switch analysis
        case 1
            behStr = 'EV';
        case 2
            behStr = 'H';
        case 3
            behStr = 'Info';
    end

    [p,h,stat] = signrank(oitM(:,1,analysis)-oitM(:,2,analysis));
    sprintf([behStr,' chosen-unchosen, EXPLOIT : p = %0.4f, sign = %2.0f'],p,stat.signedrank)
    sprintf('mean diff = %2.4f; range = [%2.4f, %2.4f]',mean(oitM(:,1,analysis)-oitM(:,2,analysis)),...
        min(oitM(:,1,analysis)-oitM(:,2,analysis)),max(oitM(:,1,analysis)-oitM(:,2,analysis)))

    [p,h,stat] = signrank(oreM(:,1,analysis)-oreM(:,2,analysis));
    sprintf([behStr,' chosen-unchosen, EXPLORE : p = %0.4f, sign = %2.0f'],p,stat.signedrank)
    sprintf('mean diff = %2.4f; range = [%2.4f, %2.4f]',mean(oreM(:,1,analysis)-oreM(:,2,analysis)),...
        min(oreM(:,1,analysis)-oreM(:,2,analysis)),max(oreM(:,1,analysis)-oreM(:,2,analysis)))

    subplot(2,3,analysis);
    title([behStr,', rule'])
    subplot(2,3,analysis+3);
    title([behStr,', explore'])
end

% compare magnitude of info and reward in exploration
sprintf('percent change in EV (chosen/unchosen) = %0.2f; in Info = %0.2f',...
    nanmean((oreM(:,1,1)-oreM(:,2,1))./oreM(:,2,1)), nanmean((oreM(:,1,3)-oreM(:,2,3))./oreM(:,2,3)))

% we can plot these as time series if we want to
% blockChange_magnitudeEnergy(data,'chosenEV',false,'newImages');
% blockChange_magnitudeEnergy(data,'chosenH',false,'newImages');
% blockChange_magnitudeEnergy(data,'chosenInfo',false,'newImages');
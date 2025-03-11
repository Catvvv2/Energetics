function dissociateOnsetExploreAndReward_magnitudeEnergy(data,behOI,interaction)

    if nargin < 3
        interaction = false;
    end

    % so reward and exploration *both* predict everything
    % we want to know which of these is the best predictor, if there is an
    % interaction between them, etc;
    y = [data.(behOI)]';
    % y = [data.magnitude_change]'; % reward influences magnitude change,
    % explore does not; explore is related to absolute magnitude though
    
    y = zscore(y); % convert to units of STD

    ores = [data.explore]';
    onsets = [data.onsets]';
    ores(onsets==1) = deal(0);
    
    X = [vertcat(data.last_rwd), ores, onsets];
    if interaction
        X = [X, X(:,1).*X(:,2),X(:,1).*X(:,3)]; % add an interaction, this is non-significant
    end

    tmp = arrayfun(@(f) repmat(f,1,length([data(f).(behOI)])),[1:length(data)],'UniformOutput',false);
    tmp = dummyvar([tmp{:}]);
    X = [X, tmp(:,1:end-1)];
    
    [b,t,stat] = glmfit(X,y);

    sprintf([behOI,'\n last rwd, main effect = %2.2f, p = %2.4f \n other ore, main effect = %2.2f, p = %2.4f \n onset, main effect = %2.2f, p = %2.4f'],...
        stat.beta(2),stat.p(2),stat.beta(3),stat.p(3),stat.beta(4),stat.p(4))

    if interaction
        sprintf(['oreXrwd interaction = %2.2f, p = %2.4f \n onsetXrwd interaction = %2.2f, p = %2.4f'],...
            stat.beta(5),stat.p(5),stat.beta(6),stat.p(6))
    end
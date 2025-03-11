function dissociateStateFromReward_magnitudeEnergy(data,behOI,predOI)

    if nargin < 3
        predOI = 'explore';
    end

    % so reward and exploration *both* predict everything
    % we want to know which of these is the best predictor, if there is an
    % interaction between them, etc;
    y = [data.(behOI)]';
    % y = [data.magnitude_change]'; % reward influences magnitude change,
    % explore does not; explore is related to absolute magnitude though
    
    y = zscore(y); % convert to units of STD

    X = [vertcat(data.last_rwd), [data.(predOI)]'];
    X = [X, X(:,1).*X(:,2)]; % add an interaction, this is non-significant
    
    tmp = arrayfun(@(f) repmat(f,1,length([data(f).(behOI)])),[1:length(data)],'UniformOutput',false);
    tmp = dummyvar([tmp{:}]);
    X = [X, tmp(:,1:end-1)];
    
    [b,t,stat] = glmfit(X,y);

    sprintf([behOI,'\n last rwd, main effect = %2.2f, p = %2.4f \n ',predOI,', main effect = %2.2f, p = %2.4f \n interaction = %2.2f, p = %2.4f'],...
        stat.beta(2),stat.p(2),stat.beta(3),stat.p(3),stat.beta(4),stat.p(4))
function fH = barrierHeightByRwdHistoryFilter_humanManifolds

home = pwd;
cd('/Users/becket/Documents/MATLAB/humanManifolds')
load('manifoldDataRaw.mat');
load('manifoldData.mat','needforCognition','Age_recode')

%%

nLags = 10; % for kernel

opts = fitoptions('Method','NonlinearLeastSquares',...
    'StartPoint',[2]);
m1 = fittype('-exp(x-b)'); % 1 param exp
xpos = [-nLags:-1]';

    data(1).T

for s = 1:length(data)

    % get the reward history filters
    rwd = data(s).reward(1:end-1); % last reward
    onsets = diff([data(s).explore])==1;
    b = wienerFilter(rwd,onsets,nLags-1,'causal');
    [f1,gof1] = fit(xpos,b,m1,opts);
    if gof1.rsquare > 0
        data(s).rwdFilter = f1.b;

        % now the thermodynamic stuff
        T = data(s).T;
        try
            pi = stationaryDist(T);
            data(s).stationaryDist = [pi(1),1-pi(1)];
            [E,~,E_b] = twoWellPotential_v2([pi(1),1-pi(1)],[1-T(1,1),T(2,1)],false);
        %     [E,~,E_b] = twoWellPotential_v2([1-pi(1),pi(1)],[1-T(2,2),1-T(1,1)],false)
            data(s).barrierHeight = E_b;
            data(s).basinDifference = E(1)-E(end); % this is delta G
        catch
            data(s).barrierHeight = NaN;
            data(s).basinDifference = NaN;
        end
    else
        data(s).rwdFilter = NaN;
        data(s).barrierHeight = NaN;
        data(s).basinDifference = NaN;
    end
end

%%
[r,p] = corrcoef([data.rwdFilter],[data.barrierHeight],'rows','complete');

fH = figure();
plot([data.rwdFilter],[data.barrierHeight],'.k')
olsline;

set(gca,'FontSize',16);
xlabel('reward filter parameter')
ylabel('barrier height')

%% This is nice and maybe relevant to the age paper
% -- there's a nice relationship between barrier height and rwd filter
% and age is part of this somehow

% dependent = needforCognition;
dependent = Age_recode;

[r,p] = corrcoef(dependent,[data.barrierHeight],'rows','complete')

fH = figure();
plot([data.barrierHeight],dependent,'.k')
% olsline;

set(gca,'FontSize',16);
ylabel('dependent')
xlabel('barrier height')

%%
cd(home);

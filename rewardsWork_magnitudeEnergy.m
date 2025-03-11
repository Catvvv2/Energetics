function [fH,data] = rewardsWork_magnitudeEnergy(data)

global params

fH = figure('Position',[476   701   500   350]);
modelStr = 'binomial';

%% first, we will ask if the neural measures mediate the relationship between rewards 

for analysis = 1:2
    switch analysis
        case 1
            behOI = 'trialSpikeRate';
        case 2
            behOI = 'magnitude';
    end

    out = NaN(length(data),5);
    for k = 1:length(data)
    
        % first, pull out the behavior and normalize it
        tmpBeh = data(k).(behOI)';
        tmpBeh = zscore(tmpBeh);
       
        % identify all the onsets
        y = [data(k).onsets]'==1;
    
        % omitted rewards predicting ore transitions, Y ~ X
        [b1,t,stat1] = glmfit(1-[data(k).last_rwd],y,modelStr);
        out(k,1) = b1(2); % keep the direct effect, this is aka "c"
        
        % mediation term, rewards predicting neural data, M ~ X
        [b2,t,stat2] = glmfit(1-[data(k).last_rwd],tmpBeh);
        
        % both predicting onsets, Y ~ X + M
        [b3,t,stat3] = glmfit([1-data(k).last_rwd,tmpBeh],y,modelStr);
        
        % now calculate the effect of the mediation
        b = b3(3);
        b_se = stat3.se(3); % mediator to dependent
        a = b2(2);
        a_se = stat2.se(2); % independent to mediator
        
        beta_indirect = a*b; % product of those two effects
        out(k,2) = beta_indirect;
        
        % can also calculate c (direct), c' (remaining direct - indirect path)
        c_prime = b3(2);
        out(k,3) = c_prime;
        
        % Goodman 1960 unbiased estimate:
        % sigma_ab = sqrt(a^2*b_se^2+b^2*a_se^2-a_se^2*b_se^2);
        
        % Sobel 1982 first order
        sigma_ab = sqrt(a^2*b_se^2+b^2*a_se^2);
        
        zprime = (a*b) ./ sigma_ab;
        % this actually is pretty good, given this distribution is nonnormal
        
        out(k,4) = 1-tcdf(zprime,length(y)-1); % save p value

        % finally, we need to check to make sure that a, b, and c are sig.
        out(k,5) = and(and(stat1.p(2) < 0.05, stat2.p(2) < 0.05), stat3.p(3) < 0.05);
    
    end
    
    % now plot
    subplot(2,2,analysis); hold on;
    for k = 1:length(data)
        bar(k,out(k,2),'LineStyle','none','FaceColor',params.colorMap(k,:));
    end
    % h(2) = bar(out(:,[3]));
    
    % test for significant mediation = sig paths + sig mediation
    tmp = find(and(out(:,4) < 0.05,out(:,5) == 1));
    plot(tmp,ones(size(tmp)).*1.4,'*k','MarkerSize',2)
    tmp = find(and(out(:,4) < 0.005,out(:,5) == 1));
    plot(tmp-0.2,ones(size(tmp)).*1.4,'*k','MarkerSize',2)
    plot(tmp+0.2,ones(size(tmp)).*1.4,'*k','MarkerSize',2)
    
    set(gca,'FontSize',14,'XTick',[1:length(data)])
    ylabel({behOI,'mediation effect'})
    xlabel('session #')
    ylim([-1.5 1.5])
    
end

%% then we will do a Wiener filter analysis

subplot(2,2,3); hold all;

xpos = [-params.nLagsWiener:-1]';
xeval = [-params.nLagsWiener:0.1:-1]';
keepB = [];

if params.nParamsWiener == 1
    % 1 parameter exponential function
    opts = fitoptions('Method','NonlinearLeastSquares',...
        'StartPoint',[2]);
    m1 = fittype('-exp(x-b)');
else
    opts = fitoptions('Method','NonlinearLeastSquares',...
        'StartPoint',[0,2]);
    m1 = fittype('a-exp(x-b)'); % proper 2 param exp% proper 2 param exp
end

for f = 1:length(data)
    b = wienerFilter(data(f).last_rwd(2:end),data(f).onsets(2:end),params.nLagsWiener-1,'causal');
    plot(xpos,b,'.','Color',params.colorMap(f,:),'MarkerSize',10);

    % now curve fitting:
    [f1,gof1] = fit(xpos,b,m1,opts);
    plot(xeval,f1(xeval),'-','Color',params.colorMap(f,:),...
        'LineWidth',1);

    if isfield(f1,'a')
        data(f).rwdFilterA = f1.a;
    end
    data(f).rwdFilterB = f1.b;
    keepB = [keepB,b];
end

b = nanmean(keepB,2);
plot(xpos,b,'.','Color','k','MarkerSize',25);
[f1,gof1] = fit(xpos,b,m1,opts);
plot(xeval,f1(xeval),'-','Color','k',...
    'LineWidth',3);

xlim([-(params.nLagsWiener+1) 0])
set(gca,'FontSize',14)
ylabel('outcome weight')
xlabel('previous trial')

% then plot the filters
subplot(2,2,4); hold all;
for k = 1:length(data)
    bar(k,data(k).rwdFilterB,...
        'LineStyle','none','FaceColor',params.colorMap(k,:));
end
ylim([0,3])
set(gca,'XTick',[1:length(data)],'FontSize',14)
xlabel('session #'); ylabel('rwd filter parameter')

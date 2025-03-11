function fH = baselinePlots_magnitudeEnergy(data)

global params
params.exSession = 3;
params.exSequence = 449:798;

fH = figure(99);
set(fH,'Position',[476   552  854   556]);

xpos = 1:length(params.exSequence);

imChoices = 1-data(params.exSession).imChoice(params.exSequence);
locChoices = 1-data(params.exSession).locChoice(params.exSequence);
best = data(params.exSession).choseBest(params.exSequence);

subplot(2,7,1:5); hold all;
plot(xpos,imChoices.*.6+1.2,'.-');
plot(xpos,locChoices.*.6+0.2,'.-');

% now illustrate the blocktype:
blocks = cumsum(data(params.exSession).blockChange(params.exSequence));
blockTypes = data(params.exSession).blockType(params.exSequence);
wonk = [];
for bl = 0:max(blocks)
    type = blockTypes(blocks == bl);
    if type(1) == 1 % what
        sType = imChoices(and(blocks == bl,best));
    elseif type(1) == 2 % where
        sType = locChoices(and(blocks == bl,best));
    end
    % now concatenate and re-code so it better matches the state labels
    wonk = [wonk; [2-type+1, repmat(sType(1),size(type))]];
end

h = plot(xpos,(wonk(:,1)-1)+wonk(:,2).*0.8 + 0.1,'s');
set(h,'MarkerFaceColor',get(h,'Color'));
% legend({'image','location'},'Location','northoutside')

h = plot(xpos,-data(params.exSession).stateLabels(params.exSequence)./14 + 2.4,'s');
set(h,'MarkerFaceColor',get(h,'Color'));

set(gca,'FontSize',16,...
    'YTick',[0.5,1.5],...
    'YTickLabel',{'location','image'},...
    'YTickLabelRotation',90);
xlabel('trial number');

% 1 what, 2 where
% plot(params.exSequence,,'s');

ylim([0, 2.4])

%% now I want to know the agreement between the HMM labels and the underlying rules

u_states = 1:5;
u_rules = 1:4;

out = NaN(length(u_states),length(u_rules),length(data));

for f = 1:length(data)
    % now illustrate the blocktype:
    blocks = cumsum(data(f).blockChange);
    blockTypes = data(f).blockType;
    best = data(f).choseBest;
    imChoices = data(f).imChoice;
    locChoices = data(f).locChoice;

    wonk = [];
    for bl = 1:max(blocks)
        type = blockTypes(blocks == bl);
        if type(1) == 1 % what
            sType = imChoices(and(blocks == bl,best));
            if isempty(sType)
                sType = 1-imChoices(and(blocks == bl,~best));
            end
        elseif type(1) == 2 % where
            sType = locChoices(and(blocks == bl,best));
            if isempty(sType)
                sType = 1-locChoices(and(blocks == bl,~best));
            end
        end
        wonk = [wonk; [type, repmat(sType(1),size(type))]];
        % 1 = what, 2 = where
    end
    
    rules = (wonk(:,1) + wonk(:,2)./2).*2 - 1;
    states = data(f).stateLabels';

    [xx,yy] = meshgrid(u_rules,u_states);
    out(:,:,f) = arrayfun(@(x,y) nanmean(rules(states==y)==x),xx,yy);

end

subplot(2,6,7:8); colormap(flipud(gray));
imagesc(nanmean(out,3));
colorbar;
set(gca,'FontSize',16,...
    'XTick',[1:4],'XTickLabel',{'im1','im2','loc1','loc2'},...
    'YTick',[1:5],'YTickLabel',{'ore','im1','im2','loc1','loc2'})

xlabel('rewarded rule')
ylabel('HMM state')


nanmean(out,3)
nanstd(out,[],3)./sqrt(size(out,3)-1)

% summarize some stuff
% first, the diagonal == how likely to match
tmp = arrayfun(@(k) diag(out(2:end,:,k)),1:size(out,3),'UniformOutput',false)
mean(mean([tmp{:}]))
[min(mean([tmp{:}])) max(mean([tmp{:}]))]

intra = [0, 1, 0, 0; 1, 0, 0, 0; 0, 0, 0, 1; 0, 0, 1, 0];
extra = [0, 0, 1, 1; 0, 0, 1, 1; 1, 1, 0, 0; 1, 1, 0, 0];
intra = arrayfun(@(k) sum(sum(intra.*out(2:end,:,k)))./sum(sum(intra)), 1:size(out,3));
extra = arrayfun(@(k) sum(sum(extra.*out(2:end,:,k)))./sum(sum(extra)), 1:size(out,3));
m = [intra; extra];
nanmean(extra-intra)
[min(extra-intra),max(extra-intra)]

% figure('Position',[476   642   219   224]); hold all;
subplot(2,9,14:15); hold all;
for f = 1:length(data)
    h = plot([0,1],m(:,f),'.-','MarkerSize',20);
    set(h,'Color',params.colorMap(f,:));
end

set(gca,'FontSize',params.FontSize,...
    'XTick',[0 1],'XTickLabel',{'intra','extra'})
xlabel('\Delta dimension');
ylabel('p(confused)')
xlim([-0.5 1.5]);

[p,h,stat] = signrank(intra,extra)

%%

theta = NaN(length(data),3);
lik = NaN(length(data),2);
times = []; ns = [];

for f = 1:length(data)
    l_times = [];

    % image choice, then location choice
    switches = diff(data(f).choices)~=0;

    for col = 1:size(switches,2)
        l_times = [l_times; diff(find(switches(:,col)))];
    end

    [theta(f,:),~,lik(f,1)] = exp2mix(l_times-1);
    times = [times; l_times];
    ns = [ns; length(l_times)];

    % now a constrained model that keeps the short at 1.5
    % and requires that there be exactly equal switching
    [~,~,lik(f,2)] = exp2mix_constrained(l_times-1);
end

subplot(2,10,18:20);
twoMixLogPlot(times-1)

% now for comparison w/ the restricted model
[h,p,stat] = lratiotest(lik(:,1),lik(:,2),2) % we can do a lik ratio test
% dof is the number of restrictions: so we've eliminated 2 parameters

subplot(2,10,9:10); hold on;
for f = 1:length(lik)
    h = plot([0,1],lik(f,:),'.-','MarkerSize',20);
    set(h,'Color',params.colorMap(f,:));
end
set(gca,'FontSize',16,'XTick',[0,1],...
    'XTickLabel',{'4 state','5 state'})
xlim([-0.5, 1.5])
ylabel('-log lik')


% if there is only random/fast switching in 1 dimension when there is
% organization/repetition in the other, then the half life should be:
%        1 + 1/2 = 1.5    (compare to theta(1))
% and we should be able to calculate the frequency from the long half life
% exactly as:

% % expected % of time in each state???
% ((theta(:,2)+1)./(theta(:,1)+1)) ./ ((theta(:,1)+1) + (theta(:,2)+1))

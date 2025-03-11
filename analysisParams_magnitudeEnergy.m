%% analysis params, magnitude + energy project
    % to be used in combination with analyze_magnitudeEnergy

global params % declare global

% which file to analyse?
params.dataFile = 'averbeckWhatWhere.mat';
params.figFormat = '-dpng'; % '-depsc'
params.colorMap = [1.0000    0.4364    0.5348;...
    0.9041    0.4925    0.3099;...
    0.6884    0.5697    0.1799;...
    0.3055    0.6804    0.2114;...
    0.1799    0.6824    0.5619;...
    0.1985    0.6419    0.9077;...
    0.4147    0.5397    1.0000;...
    0.8634    0.4285    1.0000];

% for block-change aligned plots, i.e. blockChange_magnitudeEnergy
params.minTrial = -5; % also used as the baseline period for sig. testing
params.maxTrial = 15;

% get the data directory we want to use
params.dataDir = '/Users/becket/Documents/MATLAB/magnitudeEnergy/';
params.figDir = '/Users/becket/Documents/MATLAB/magnitudeEnergy/figures';
addpath(params.dataDir);

% plotting parameters
params.FontSize = 16;
params.defaultBlack = [0.1 0.1 0.1]; % not quite black
params.plotRaw = true; % generic plot raw logical
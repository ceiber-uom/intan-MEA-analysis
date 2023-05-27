
clear
% list = dir('../*/*.rhd');

ks_path = ['\\shared.sydney.edu.au\research-data\PRJ-vnrg2019\' ...
                              'V19_transfer\Calvin\MEA\TestKS\'];

list = dir([ks_path '/*.rhd']);

p_ = @(x) [x.folder filesep x.name]; % path expander

data = tools.readIntan(p_(list(1)));

%%


data = tools.spikeDetection(data); 


%%

plots.epochs(data);


%% Demonstrate the power of the line noise filter

cla, hold on
plot(data.AMP.time, data.AMP.data(:,24));

data.AMP.data = tools.removeLineNoise(data.AMP.data', ...
                            data.config.general.amplifier_sample_rate)';

plot(data.AMP.time, data.AMP.data(:,24));

%% 

data = tools.removeCommonMode(data);
data = tools.segmentEpochs(data);

%%

plots.epochs(data);

%%


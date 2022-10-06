
clear
list = dir('../05*/*.rhd');

p_ = @(x) [x.folder filesep x.name]; % path expander

data = tools.readIntan(p_(list(1)));

%% Demonstrate the power of the line noise filter

cla, hold on
plot(data.AMP.time, data.AMP.data(:,24));

data.AMP.data = tools.removeLineNoise(data.AMP.data', ...
                            data.config.general.amplifier_sample_rate)';

plot(data.AMP.time, data.AMP.data(:,24));

%% 

data = tools.segmentEpochs(data);

%%


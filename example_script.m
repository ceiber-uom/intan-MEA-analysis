
clear
list = dir('../05*/*.rhd');

p_ = @(x) [x.folder filesep x.name]; % path expander

data = tools.readIntan(p_(list(1)));

%%
data = tools.segmentEpochs(data);


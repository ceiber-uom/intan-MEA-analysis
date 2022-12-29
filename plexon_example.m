clear
list = dir('../Pro*/*.xlsx'); % .csv
p_ = @(x) [x.folder filesep x.name]; % path expander

data = tools.readPlexon(p_(list(1)));

% tw = load(strrep(p_(list(1)), '.csv','.mat'));

%%

msa = analysis.multiState(data, '-roi', [0 500], ...
                                '-chan',[56 57 60 62], ... 
                                '-bin', 0.5);

%% 

all_ba = analysis.burstAnalysis(data); 

% ba = analysis.burstAnalysis(data, '-chan', 57, '-unit', 4);
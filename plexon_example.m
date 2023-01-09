

clear

list = dir('../Pro*/*.xlsx'); % .csv
p_ = @(x) [x.folder filesep x.name]; % path expander

data = tools.readPlexon(p_(list(1)));

% tw = load(strrep(p_(list(1)), '.csv','.mat'));

%%

msa = analysis.states(data, '-roi', [0 500], ...
                            '-chan', 33:64,  ... 
                            '-bin', 0.5);

%% 

ba = analysis.bursts(data); 

% ba = analysis.burstAnalysis(data, '-chan', 57, '-unit', 4);

%%

plots.STA(data,'-trigger',[38 1],'-label')

clear

ks_folder = ['\\shared.sydney.edu.au\research-data\PRJ-vnrg2019\' ...
                            'V19_transfer\Calvin\MEA\Fl_int_durNF\'];

%%

list = dir([ks_folder '/*.rhd']); % list of multiple files
data = tools.readIntan(list,'ADC',1); % read multiple files and merge

%%
data = tools.readKS(data,'-quality',4); % read ks_sorted spike file


%%

ed = tools.segmentEpochs(data); % epoched data

%%

plots.raster(ed)
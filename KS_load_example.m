
clear

ks_folder = ['\\shared.sydney.edu.au\research-data\PRJ-vnrg2019\' ...
                            'V19_transfer\Calvin\MEA\Fl_int_durNF\'];

list = dir([ks_folder '/*.rhd']); % list of multiple files
data = tools.readIntan(list,'ADC',1); % read multiple files and merge
data = tools.readKS(data,'-quality',4); % read ks_sorted spike file

%%
ed = tools.segmentEpochs(data); % epoched data

%%

% plots.raster(ed); % ,'-per-unit'

clf
plots.psth(ed) % ,'-per-unit' to see every unit on a seperate axis


%% Zoom into a set of units on a pair of channels (arbitrarially chosen)

clf
plots.psth(ed,'-unit',90:100,'-label','-chan',13:14)
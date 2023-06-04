
clear

ks_folder = ['\\shared.sydney.edu.au\research-data\PRJ-vnrg2019\' ...
                            'V19_transfer\Calvin\MEA\Fl_int_durNF\'];

list = dir([ks_folder '/*.rhd']); % list of multiple files
data = tools.readIntan(list,'ADC',1); % read multiple files and merge
data = tools.readKS(data,'-quality',4); % read ks_sorted spike file

epoched = tools.segmentEpochs(data); % epoched data

%%

% plots.raster(ed); % ,'-per-unit'

clf
% plots.psth(ed)

plots.psth(epoched,'-per-unit'); %  to see every unit on a seperate axis


%% Zoom into a set of units on a pair of channels 

% get(gco,'userdata') % return the [channel unit_id pass_id] of a selected histogram

clf
plots.psth(epoched,'-per-unit','-unit',17:22,'-chan',16:21); %  see just these units


%%

plots.response_curve(epoched,'-roi',[0 3]);

%%

plots.psth(epoched,'-roi', [-0.1 2.5]); %  to see every unit on a seperate axis


% todo - implement epochs.block_id for averaging 
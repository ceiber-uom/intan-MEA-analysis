
clear

% ks_folder = ['\\shared.sydney.edu.au\research-data\PRJ-vnrg2019\' ...
%                             'V19_transfer\Calvin\MEA\Fl_int_durNF\'];
% list = dir([ks_folder '/*.rhd']); % list of multiple files

data = tools.readIntan('','ADC',1); % read one or multiple files and merge
data = tools.readKS(data,'-quality',4); % read ks_sorted spike file

% data = tools.segmentEpochs(data); % epoched data
% setting data.epochs.condition_id enables averaging in the PSTH 

epochs = tools.segmentEpochs(data,'-jitter',0.1); % epoched data
epochs = set_condition_id(epochs); % for 01_ stimuli

% data.epochs.condition_id = [ones(50,1); 2*ones(50,1)]; % 02_SQ1X100NF.rhd
% data.epochs.condition_id = ones(50,1); %  03_SQ3x50NF.rhd 

%% 

% plots.raster(epochs); % ,'-per-unit'

clf
% plots.psth(epochs) % one panel per channel
% plots.psth(data,'-per-unit'); %  to see every unit on a seperate axis

% 03_SQ3x50NF.rhd - zoom into a relevent ROI for this stimulus:
plots.psth(epochs,'-per-unit','-time',1,'-roi',[-0.01 0.05],'-sem');

% try adding -sem for an alternate plot style


%% Zoom into a set of units on a pair of channels 

% get(gco,'userdata') % return the [channel unit_id pass_id] of a selected histogram

clf
plots.psth(data,'-per-unit','-unit',17:22,'-chan',16:21); %  see just these units

% the -chan is necessary because unit_ids are not assumed to be unique
% across channels (i.e. the way BlackRock or Plexon handles unit_its). i.e.
% Ch1.u2 and Ch10.u2 are not the same cell. This is different to KiloSort
% which assigns unique unit_ids across cells (which might be picked up on
% multiple channels), so Ch1 might have units 1-4, Ch2 4-6, Ch3 5-12, etc. 


%%

plots.response_curve(epochs,'-x','duration','-roi',[-2 5],'-per-unit');

%%

plots.psth(data,'-roi', [-0.1 2.5]); %  to see every unit on a seperate axis


% todo - implement epochs.block_id for averaging 
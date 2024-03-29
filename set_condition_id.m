
function stim = set_condition_id(stim, protocol) %#ok<INUSD> 
% this was developed quickly and is a special case file at the moment.
% Normally this is an example of something that would be 'user script'
% code as opposed to 'engine code'. 
% 
% In theory this can work for any amplitude-duration grid but might need
% to be tweaked for different files. 

% assert protocol = 01_SQ1100NF565040423

if isfield(stim,'epochs') % enable call on top-level object
    stim.epochs = set_condition_id(stim.epochs);
    return
end

% you can see if this is 'correct' as stim.duration and stim.average should
% form a grid (potentially with non-uniform spacing):
% >> clf, plot(stim.duration, stim.average, 'o')

stim.duration = round(stim.duration,3); % correct noise in estimate
stim.average  = round(stim.average,2);  % correct noise in estimate

% for this protocol, 
% just rounding works for all but the lowest amplitude stimulus pulse
a_levels = unique(stim.average(stim.duration == max(stim.duration)));
fix_me = (stim.duration == min(stim.duration)); % these need more help

a_noisy = stim.average(fix_me); % .average for shortest pulses
a_noisy(a_noisy <= min(a_levels)) = min(a_levels); 

for ii = 2:numel(a_levels) % round up to next level
    a_noisy(a_noisy <= a_levels(ii) & ...
            a_noisy > a_levels(ii-1)) = a_levels(ii);
end

stim.average(fix_me) = a_noisy; % put these back in the stim data structure

[~,~,cond_id] = unique([stim.duration stim.average],'rows');
n_reps = arrayfun(@(c) sum(cond_id == c), 1:max(cond_id)); 

assert(all(n_reps > 1),'please double-check stimulus counts')

stim.condition_id = cond_id; 

return
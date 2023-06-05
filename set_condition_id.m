

function stim = set_condition_id(stim, protocol) %#ok<INUSD> 
% assert protocol = 01_SQ1100NF565040423

stim.duration = round(stim.duration,3); % correct noise in estimate
stim.average  = round(stim.average,2);  % correct noise in estimate

% just rounding works for all but the lowest amplitude stimulus pulse

a_levels = unique(stim.average(stim.duration == max(stim.duration)));

fix_me = (stim.duration == min(stim.duration)); % these need more help

a_noisy = stim.average(fix_me); 
a_noisy(a_noisy <= min(a_levels)) = min(a_levels); 

for ii = 2:numel(a_levels)

    a_noisy(a_noisy <= a_levels(ii) & ...
            a_noisy > a_levels(ii-1)) = a_levels(ii);
end

stim.average(fix_me) = a_noisy;

[~,~,cond_id] = unique([stim.duration stim.average],'rows');
n_reps = arrayfun(@(c) sum(cond_id == c), 1:max(cond_id)); 

assert(all(n_reps > 1),'please double-check stimulus counts')

stim.condition_id = cond_id; 

return


function spikes = simplify(data, varargin)
% [spikes] = tools.simplify( data, ... ) 
% 
% convert the format of the supplied data.SPIKE from a 'flattened'
% structure to a more user-friendly format. The output struct array has one
% entry per channel / unit, and is formatted as follows: 
% 
% spikes(x).time    - cell array of spike-times, one per pass (see -flat) 
% spikes(x).channel - origin spike channel_id
% spikes(x).unit    - origin spike unit_id
% spikes(x).shape   - shape variables for each spike (as exported by
%                      Plexon, see data.config.shape_variables)
% 
% ... other fields for epoched data ? 
%
% Options:
%  -flat  Flatten spikes to one vector across all passes (for epoched data) 
%          (default: output a cell array of spike-times for each pass)
%  -undo  Reverse this data transformation on the input struct array
% 
% (inherited from tools.forChannels):
%  -chan [c1 c2 c3 ...] Set channels 
%  -unit [u1 u2 ... ]   Filter for the specified units 
%  -roi  [t0 t1]        Filter spikes in time
%  -merge-units         Ignore unit codes
%  -pass [pass_ids]     Select passes to keep (for expoched data)
%  -no-hash             Only show spikes with unit code > 0. 
%
% for Dario
% v0.1 - Calvin Eiber <c.eiber@ieee.org> 1 July 2023


named = @(s) strncmpi(s,varargin,numel(s));

do_pass_vector = any(named('-flat')); 

if any(named('-undo')), spikes = undo_unpacking(data); return, end

spk_data = @(d,varargin) (d); % extract the data

printInfo(); 

[~, spikes] = tools.forChannels(data, spk_data, varargin{:},'--ordered'); 

if isfield(spikes,'pass_begin')
     n_passes = numel(spikes(1).pass_begin);
     t0 = [0; spikes(1).pass_begin];

     error needs_testing

else n_passes = max(cat(1,spikes.pass));
     t0 = zeros(n_passes+1,1);
end

spk_times = @(p,s) reshape(s.time(s.pass == p),1,[])-t0(s.pass+1); 

for ii = 1:numel(spikes)

    spikes(ii).channel = spikes(ii).channel(1);
    spikes(ii).unit = spikes(ii).unit(1);

    if do_pass_vector && n_passes > 1
        get_st = @(p) spk_times(p,spikes(ii)); 
        spikes(ii).time = arrayfun(get_st, 1:n_passes,'unif',0);
    end
end
disp('Done! ')

if ~do_pass_vector, spikes = rmfield(spikes,'pass'); end

%%
if nargout == 0
    assignin('caller','spikes',spikes), clear
end
return




function d = undo_unpacking(spk)
% undo function for the above unpacking. 
% put the spikes back into a single structure. 
% used a lot by my plot code which expects single structure spikes

fix_passes = iscell(spk(1).time);

if fix_passes, 
    n_pass = numel(spk(1).time); 
    [spk.pass] = deal([]); 
    pass_ids = @(p,s) p*ones(size(s.time{p}));

end

for ii = 1:numel(spk)

    if fix_passes,
        spk(ii).pass = arrayfun(@(p) pass_ids(p,spk(ii)), 1:n_pass);
        spk(ii).time = reshape([spk(ii).time{:}],[],1); 
    end
    spk(ii).channel = spk(ii).channel(1) * ones(size(spk(ii).time));
    spk(ii).unit    = spk(ii).unit(1)    * ones(size(spk(ii).time));
end

d = struct;
for var = {'time','channel','unit','shape'}
    if ~isfield(spk,var{1}), continue, end
    d.(var{1}) = cat(1,spk.(var{1})); 
end


return
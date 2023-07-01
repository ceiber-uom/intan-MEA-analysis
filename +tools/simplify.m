

function spikes = simplify(data, varargin)
% [spikes] = tools.simplify( data, ... ) 
% 
% convert the format of the supplied SPIKE data from a single vector with
% all spikes 
%
%
% for Dario


named = @(s) strncmpi(s,varargin,numel(s));

do_pass_vector = any(named('-pass-ok')); 

spk_data = @(d,varargin) (d); % extract the data

printInfo(); 
[~, spikes] = tools.forChannels(data.SPIKE, spk_data, varargin{:},'--ordered'); 

n_passes = max(cat(1,spikes.pass))

for ii = 1:numel(spikes)

    spikes(ii).channel = spikes(ii).channel(1);
    spikes(ii).unit = spikes(ii).unit(1);

    if do_pass_vector && n_passes > 1

        spk_times = @(p) reshape(spikes(ii).time(spikes(ii).pass == p), ...
                                1,[]); % - pass(p).time ? 
        spikes(ii).time = arrayfun(spk_times, 1:n_passes,'unif',0);
    end
end
disp('Done! ')

if ~do_pass_vector, spikes = rmfield(spikes,'pass'); end

%%
if nargout == 0
    assignin('caller','spikes',spikes), clear
end
return
function results = f012(data, varargin)
% results = analysis.f012( spike_data, ... )
% 
% For each channel and unit in the supplied data, this analysis looks at
%  the stimulus-response relationship. This requires epoched data. The
%  following data structure is returned for each channel / unit : 
% 
%  .channel_unit   : [channel unit] of analysis
%  
%  <<< to be implemented >>> 
% 
% Options: 
%  -pdf              : generate a PDF documenting the analysis outcomes
%  -roi [start end]  : set time analysis window (default: whole wave)
% 
% (inherited from tools.forChannels)
% -chan [c1 c2 ... ] : Select channels to analyse
% -unit [u1 u2 ... ] : Filter for the specified units 
% -merge-units       : Combined analysis of all spikes on each channel
% -pass [pass_ids]     Filter for epoch IDs to analyse
% -no-hash             Only analyse spikes with unit code > 0. 
% 
% v0.1 - 28 December 2022 - Calvin Eiber <c.eiber@ieee.org>

if nargin == 0, try data = evalin('caller','data'); end, end %#ok<TRYNC> 

disp(datestr(now)), disp('Running f012 analysis')

do_PDF = any(named('-pdf'));
plots.PDF_tools('setup',do_PDF);

printInfo(); 
[~, results] = tools.forChannels(data, @run_f012, varargin{:},'--ordered'); 
disp('Done! ')

% make_summary_graphic(results, varargin)

plots.PDF_tools('compile',do_PDF,'f012-analysis (%d).pdf')

return

%% per-channel / per-unit analysis
function stats = run_f012(data, index, varargin)

printInfo('f012 analysis (c%d.u%d) ', index);

named = @(s) strncmpi(s,varargin,numel(s));
get_ = @(v) varargin{find(named(v))+1};

error implement_f012_analysis


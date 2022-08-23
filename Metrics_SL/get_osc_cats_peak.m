function [cats] = get_osc_cats_peak(NumPeaks,OffTimes,varargin )
% Calculates the percentage of trajectories that fall into different
% oscillation categories in the time doc using number of peaks as criterium
%----------------------------------------------------------------------
% [cats, frac] = get_osc_cats(NumPeaks,OffTimes,varargin )
%--------------------------------------------------------------------------
% INPUT:
%   REQUIRED:
%       NumPeaks:           number of peaks
%       OffTimes:          time of signal termination
%   OPTIONAL:
%       'Threshold':          cutoff for oscillatory & non-oscillatory
%       (default: 3)

% OUTPUT:
%   tbl:     fraction of cells that are 'on', 'oscillatory', 'transient',
%           'persistent'

p=inputParser; 
addRequired(p, 'NumPeaks', @isnumeric); 
addRequired(p, 'OffTimes', @isnumeric); 
addParameter(p, 'Threshold', 3, @isnumeric); 
parse(p,NumPeaks,OffTimes,varargin{:}); 
NumPeaks = p.Results.NumPeaks; 
OffTimes=p.Results.OffTimes;
cats = repmat({'non_osc'}, size(NumPeaks)); 
cats(NumPeaks>=p.Results.Threshold)={'osc'};
cats(OffTimes==0)={'off'};
names = {'off', 'osc', 'non_osc'};
cats= categorical(cats,names); 

%{
var_names = [names,'on'];
counts =array2table(countcats(cats)', 'VariableNames', categories(cats));
% counts.on = numel(cats) -counts.off;
counts.on = sum(~(OffTimes==0));
% var_names = {'on', 'osc', 'non_osc'};
% frac = counts(:,var_names);
frac = counts(:,2:end);
frac.on = round(counts.on/numel(cats),2);
frac{:,var_names(3:4)} = round(counts{:,var_names(3:4)}./counts.on,2);
%}

end
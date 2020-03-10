function [] = see_NFkB_HM_sorted_test(id,varargin)


%% INPUT PARSING
% Create input parser object, add required params from function input
p = inputParser;
% Required: ID input
valid_id = @(x) assert((isnumeric(x)&&length(x)==1)||exist(x,'file'),...
    'ID input must be spreadsheet ID or full file path');
addRequired(p,'id',valid_id);
parse(p,id, varargin{:})
%%
%[graph, info, measure] = see_nfkb_native(id, varargin)
%[graph, info, measure] = see_nfkb_native(id);
%[metrics,aux, graph, info, measure] = nfkbmetrics(id,varargin);
[metrics,aux, graph, info, measure] = nfkbmetrics(id);

[~,graph.order] = sort(metrics.peakfreq,'descend'); 
figs.a = figure('name','ColormapStack');
set(figs.a,'Position', [500 700 1200 600])
colormapStack(graph.var(graph.order,:),graph.celldata(graph.order,:), graph.opt); 
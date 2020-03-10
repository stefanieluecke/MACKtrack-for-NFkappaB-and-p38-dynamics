function [] = see_NFkB_HM_sorted_test(id, varargin)


%% INPUT PARSING
% Create input parser object, add required params from function input
p = inputParser;
% Required: ID input
valid_id = @(x) assert((isnumeric(x)&&length(x)==1)||exist(x,'file'),...
    'ID input must be spreadsheet ID or full file path');
addRequired(p,'id',valid_id);

% Optional parameters
expectedFlags = {'on','off'};
addParameter(p,'Verbose','off', @(x) any(validatestring(x,expectedFlags)));%checks whether optional name-value argument matches on or off %checks if x matches expectedFlags
valid_conv = @(x) assert(isnumeric(x)&&(x>=0)&&(length(x)==1),...
    'Parameter must be single integer >= 0'); %checks whether parameters below are single integers
addParameter(p,'ConvectionShift',1, valid_conv); %allows adjustment of convection shift (?)
addParameter(p,'MinLifetime',100, @isnumeric); %allows adjustment of minimum lifetime (?)
addParameter(p,'MinSize',90, valid_conv); %allows adjustment of minimum size (?)
addParameter(p,'StartThresh',2, valid_conv); %max allowable starting threshhold to filter out cells with pre-activated NFkB, default is 2
addParameter (p, 'Baseline', 1, @isnumeric);
addParameter (p, 'GraphLimits',[-0.25 10],@isnumeric);
addParameter(p,'TrimFrame',157, @isnumeric);
%? add Parameter check to SortMetric!
addParameter(p, 'SortMetric', 'peakfreq');
expectedOrder = {'ascend', 'descend'};
addParameter(p, 'SortOrder', 'descend', @(x)any(validatestring(x, expectedOrder))); 

parse(p,id, varargin{:})
%%

[metrics,aux, graph, info, measure] = nfkbmetrics_SL(id, 'MinLifetime',p.Results.MinLifetime,...
                            'ConvectionShift',p.Results.ConvectionShift, 'Baseline',p.Results.Baseline,...
                            'MinSize', p.Results.MinSize,'StartThresh', p.Results.StartThresh, 'Verbose', ... 
                            p.Results.Verbose, 'GraphLimits', p.Results.GraphLimits, 'TrimFrame', p.Results.TrimFrame);
SortMetric = p.Results.SortMetric;
[~,graph.order] = sort(metrics.(SortMetric), p.Results.SortOrder);
figs.a = figure('name','ColormapStack');
set(figs.a,'Position', [500 700 1200 600])
colormapStack(graph.var(graph.order,:),graph.celldata(graph.order,:), graph.opt); 
% a = metrics;
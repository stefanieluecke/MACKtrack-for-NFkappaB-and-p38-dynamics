function [] = see_NFkB_KTR_HM_sorted_test(id, varargin)


%% INPUT PARSING
% Create input parser object, add required params from function input
p = inputParser;
% Required: ID input
valid_id = @(x) assert((isnumeric(x)&&length(x)==1)||exist(x,'file'),...
    'ID input must be spreadsheet ID or full file path');
addRequired(p,'id',valid_id);

% Optional parameters
expectedFlags = {'on','off'};
addParameter(p,'Verbose','on', @(x) any(validatestring(x,expectedFlags)));%checks whether optional name-value argument matches on or off %checks if x matches expectedFlags
valid_conv = @(x) assert(isnumeric(x)&&(x>=0)&&(length(x)==1),...
    'Parameter must be single integer >= 0'); %checks whether parameters below are single integers
addParameter(p,'ConvectionShift',1, valid_conv); %allows adjustment of convection shift (?)
addParameter(p,'MinLifetime',100, @isnumeric); %allows adjustment of minimum lifetime (?)
addParameter(p,'MinSize',90, valid_conv); %allows adjustment of minimum size (?)
addParameter(p,'TrimFrame',157, @isnumeric);
addParameter(p,'StartThreshNFkB',2, valid_conv); %max allowable starting threshhold to filter out cells with pre-activated NFkB, default is 2
addParameter (p, 'BaselineNFkB', 1, @isnumeric); %? not used in code?
addParameter (p, 'GraphLimitsNFkB',[-0.25 8],@isnumeric);
addParameter(p,'StartThreshKTR',60, valid_conv); %max allowable starting threshhold to filter out cells with pre-activated KTR, default is 0.6
addParameter (p, 'BaselineKTR', 1, @isnumeric);%? not used in code?
addParameter (p, 'GraphLimitsKTR',[0 0.35],@isnumeric);
addParameter(p, 'SortMetric', 'peakfreq_nfkb');
expectedOrder = {'ascend', 'descend'};
addParameter(p, 'SortOrder', 'descend', @(x)any(validatestring(x, expectedOrder))); 

parse(p,id, varargin{:})
%%
[metrics,aux, graph, info, measure] = nfkb_ktr_metrics(id, 'MinLifetime',p.Results.MinLifetime,...
                            'ConvectionShift',p.Results.ConvectionShift, 'BaselineNFkB',p.Results.BaselineNFkB,'BaselineKTR',p.Results.BaselineKTR,...
                            'MinSize', p.Results.MinSize,'StartThreshNFkB', p.Results.StartThreshNFkB,'StartThreshKTR', p.Results.StartThreshKTR, 'Verbose', ... 
                            p.Results.Verbose, 'GraphLimitsNFkB', p.Results.GraphLimitsNFkB,'GraphLimitsKTR', p.Results.GraphLimitsKTR, 'TrimFrame', p.Results.TrimFrame);

SortMetric = p.Results.SortMetric;
[~,graph.order] = sort(metrics.(SortMetric), p.Results.SortOrder);

figs.a = figure('name','Heatmap_NFkB');
set(figs.a,'Position', [500 7 800 600])
colormapStack(graph.var_nfkb(graph.order,:),graph.celldata(graph.order,:), graph.opt_nfkb); %calls subfunction that makes figure    

figs.b = figure('name','Heatmap_KTR');
set(figs.b,'Position', [500 7 700 600])
colormapStack(graph.var_ktr(graph.order,:),graph.celldata(graph.order,:), graph.opt_ktr); %calls subfunction that makes figure    

%add in line plot of averages and small multiple plots
end

function [ID] = see_NFkB_KTR_ratio_multiple_HMs(IDs, varargin)

% 
%Enter highest responder first
%% INPUT PARSING
% Create input parser object, add required params from function input
p = inputParser;
% Required: ID input
valid_id = @(x) assert((isnumeric(x)&&length(x)==1)||exist(x,'file'),...
    'ID input must be spreadsheet ID or full file path');
addRequired(p,'IDs');%vector containing IDs to be plotted

expectedFlags = {'on','off'};
addParameter(p,'Verbose','off', @(x) any(validatestring(x,expectedFlags)));%checks whether optional name-value argument matches on or off %checks if x matches expectedFlags
valid_conv = @(x) assert(isnumeric(x)&&(x>=0)&&(length(x)==1),...
    'Parameter must be single integer >= 0'); %checks whether parameters below are single integers
addParameter(p,'ConvectionShift',1, valid_conv); %allows adjustment of convection shift (?)
addParameter(p,'MinLifetime',109, @isnumeric); %allows adjustment of minimum lifetime (?)
addParameter(p,'MinSize',90, valid_conv); %allows adjustment of minimum size (?)
addParameter(p,'TrimFrame',157, @isnumeric);
addParameter(p,'StartThreshNFkB',14, valid_conv);%max allowable starting threshhold (before baseline deduction)to filter out cells with pre-activated NFkB
addParameter (p, 'OnThreshNFkB', 3, @isnumeric); %sigma threshold for determining responders
%addParameter (p, 'OnThreshNFkB', 8, @isnumeric); %sigma threshold for determining responders
addParameter (p, 'GraphLimitsNFkB',[-0.25 7],@isnumeric);
addParameter(p,'StartThreshKTR',0.9, valid_conv); %max allowable starting threshhold to filter out cells with pre-activated KTR, default is 0.6
addParameter (p, 'OnThreshKTR', 3, @isnumeric); %sigma threshold for determining responders
addParameter (p, 'GraphLimitsKTR',[-0.02,0.35],@isnumeric);
addParameter(p, 'SortMetric', 'peakfreq_nfkb');
expectedOrder = {'ascend', 'descend'};
addParameter(p, 'SortOrder', 'descend', @(x)any(validatestring(x, expectedOrder))); 
addParameter(p, 'SortIndex', 1) %to allow sorting by multi-column metrics
addParameter(p, 'StimulationTimePoint', 13, @isnumeric)
addParameter(p, 'FramesPerHour', 12, @isnumeric)
addParameter(p, 'NFkBBackgroundAdjustment', 'on',@(x) any(validatestring(x,expectedFlags))) %option to turn off NFkB fluorescence distribution adjustment
addParameter(p,'NFkBBaselineDeduction', 'on', @(x) any(validatestring(x,expectedFlags))) %option to turn off NFkB baseline deduction
addParameter(p,'NFkBBaselineAdjustment', 'on', @(x) any(validatestring(x,expectedFlags))) %option to turn off adjusment of NFkB trajectories with correction factor for fluorescence drop derived from Mock experiments
addParameter(p,'KTRBaselineDeduction', 'on', @(x) any(validatestring(x,expectedFlags))) %option to turn off NFkB baseline deduction
addParameter(p,'KTRBaselineAdjustment', 'on', @(x) any(validatestring(x,expectedFlags))) %option to turn off adjusment of KTR trajectories with correction factor for fluorescence drop derived from Mock experiments

expectedFilters = {'none','nfkb', 'ktr', 'both'};
addParameter(p, 'FilterResponders','none', @(x) any(validatestring(x,expectedFilters)));%filter out non-responders or not
addParameter(p, 'IncludeKTR', 'on',@(x)any(validatestring(x, expectedFlags)));

addParameter(p, 'MetricFilter', 'off',@(x)any(validatestring(x, expectedFlags)));
addParameter(p, 'FilterMetric', 'max_amplitude_4h_nfkb');
addParameter (p, 'MetricThresh', 1, @isnumeric);

addParameter (p, 'RestrictCellNumber', 'off',@(x)any(validatestring(x, expectedFlags)));
addParameter (p, 'CellNumber', []);

parse(p,IDs, varargin{:})

n = numel(IDs);
ID(n).metrics = [];
ID(n).graph = [];
ID(n).info = [];
ID(n).aux = [];
ID(n).measure = [];
%% 
for i= 1:n
[ID(i).metrics,ID(i).aux, ID(i).graph, ID(i).info, ID(i).measure] = nfkb_ktr_ratio_metrics(IDs(i), 'MinLifetime',p.Results.MinLifetime,...
                            'ConvectionShift',p.Results.ConvectionShift, 'OnThreshNFkB',p.Results.OnThreshNFkB,'OnThreshKTR',p.Results.OnThreshKTR,...
                            'MinSize', p.Results.MinSize,'StartThreshNFkB', p.Results.StartThreshNFkB,'StartThreshKTR', p.Results.StartThreshKTR, 'Verbose', ... 
                            p.Results.Verbose, 'GraphLimitsNFkB', p.Results.GraphLimitsNFkB,'GraphLimitsKTR', p.Results.GraphLimitsKTR, 'TrimFrame', p.Results.TrimFrame, 'StimulationTimePoint', p.Results.StimulationTimePoint, ...
                            'FramesPerHour', p.Results.FramesPerHour, 'NFkBBaselineDeduction', p.Results.NFkBBaselineDeduction, 'NFkBBackgroundAdjustment',p.Results.NFkBBackgroundAdjustment,...
                            'NFkBBaselineAdjustment', p.Results.NFkBBaselineAdjustment,'KTRBaselineDeduction', p.Results.KTRBaselineDeduction,'KTRBackgroundAdjustment', p.Results.KTRBackgroundAdjustment,'IncludeKTR',p.Results.IncludeKTR);
end

SortMetric = p.Results.SortMetric;
FilterMetric = p.Results.FilterMetric;
for i = 1:n
    ID(i).graph.sort_metric =   ID(i).metrics.(SortMetric);  
    ID(i).graph.filter_metric =   ID(i).metrics.(FilterMetric);  
    ID(i).graph.opt_nfkb.title = ID(i).info.name;
    ID(i).graph.opt_ktr.title = ID(i).info.name;
end 

%% Filter required graph data based on NFkB and KTR responder status
switch p.Results.FilterResponders 
    case 'nfkb'
        for i = 1:n
           ID(i).graph.var_nfkb     = ID(i).graph.var_nfkb(ID(i).metrics.responder_index_nfkb == 1,:);
           if strcmpi(p.Results.IncludeKTR,'on')
               ID(i).graph.var_ktr      = ID(i).graph.var_ktr(ID(i).metrics.responder_index_nfkb == 1,:);
           end
           ID(i).graph.celldata     = ID(i).graph.celldata(ID(i).metrics.responder_index_nfkb == 1,:);
           ID(i).graph.sort_metric  = ID(i).graph.sort_metric(ID(i).metrics.responder_index_nfkb == 1,:);
           ID(i).graph.filter_metric  = ID(i).graph.filter_metric(ID(i).metrics.responder_index_nfkb == 1,:);
           
        end
    case 'ktr'
        for i = 1:n
           ID(i).graph.var_nfkb     = ID(i).graph.var_nfkb(ID(i).metrics.responder_index_ktr == 1,:);
           ID(i).graph.var_ktr      = ID(i).graph.var_ktr(ID(i).metrics.responder_index_ktr == 1,:);
           ID(i).graph.celldata     = ID(i).graph.celldata(ID(i).metrics.responder_index_ktr == 1,:);
           ID(i).graph.sort_metric  = ID(i).graph.sort_metric(ID(i).metrics.responder_index_ktr == 1,:);
           ID(i).graph.filter_metric  = ID(i).graph.filter_metric(ID(i).metrics.responder_index_ktr == 1,:);
        end
    case 'both'
        for i = 1:n
           ID(i).graph.var_nfkb     = ID(i).graph.var_nfkb((ID(i).metrics.responder_index_nfkb == 1     & ID(i).metrics.responder_index_ktr == 1),:);
           ID(i).graph.var_ktr      = ID(i).graph.var_ktr((ID(i).metrics.responder_index_nfkb == 1      & ID(i).metrics.responder_index_ktr == 1),:);
           ID(i).graph.celldata     = ID(i).graph.celldata((ID(i).metrics.responder_index_nfkb == 1     & ID(i).metrics.responder_index_ktr == 1),:);
           ID(i).graph.sort_metric  = ID(i).graph.sort_metric((ID(i).metrics.responder_index_nfkb == 1  & ID(i).metrics.responder_index_ktr == 1),:);
           ID(i).graph.filter_metric  = ID(i).graph.filter_metric((ID(i).metrics.responder_index_nfkb == 1  & ID(i).metrics.responder_index_ktr == 1),:);
        end
end
%% Filter for a metric threshold if desired (rarely)
if strcmpi(p.Results.MetricFilter,'on')
        for i = 1:n
           ID(i).graph.var_nfkb     = ID(i).graph.var_nfkb(ID(i).graph.filter_metric >= p.Results.MetricThresh,:);
           if strcmpi(p.Results.IncludeKTR,'on')
               ID(i).graph.var_ktr      = ID(i).graph.var_ktr(ID(i).graph.filter_metric  >= p.Results.MetricThresh,:);
           end
           ID(i).graph.celldata     = ID(i).graph.celldata(ID(i).graph.filter_metric  >= p.Results.MetricThresh,:);
           ID(i).graph.sort_metric  = ID(i).graph.sort_metric(ID(i).graph.filter_metric >= p.Results.MetricThresh,:);
        end
end
 
%% Uniform random selection of cells to display if desired

if strcmpi(p.Results.RestrictCellNumber,'on')
     CellNumber = p.Results.CellNumber;
    if isempty(p.Results.CellNumber) 
        for i = 1:n
            CellNumber = [CellNumber, size(ID(i).graph.var_nfkb,1)];
        end
        CellNumber = min(CellNumber);
    end

    for i = 1:n
           index = randperm(size(ID(i).graph.var_nfkb,1),CellNumber);
           ID(i).graph.var_nfkb     = ID(i).graph.var_nfkb(index,:);
           if strcmpi(p.Results.IncludeKTR,'on')
               ID(i).graph.var_ktr      = ID(i).graph.var_ktr(index,:);
           end
           ID(i).graph.celldata     = ID(i).graph.celldata(index,:);
           ID(i).graph.sort_metric  = ID(i).graph.sort_metric(index,:);
     end
end
    
%% Sort cells by sort metric
for i= 1:n
    [~,ID(i).graph.order] = sort(ID(i).graph.sort_metric(1:end, p.Results.SortIndex), p.Results.SortOrder, 'MissingPlacement', 'last');
end

%% Plot NFkB and KTR heatmaps below each other

if strcmpi(p.Results.IncludeKTR,'on')
    figs.a = figure('name', 'HMs NFkB+KTR');

    set(figs.a,'Position', [500 400 1000 600])
    fig_handle = gcf;
    tiledlayout(2,n)

    for i= 1:n
        axes.ax(i) = nexttile;
        if i < n
            ID(i).graph.opt_nfkb.Name = '';
            SortMetricName = '';
        else
            SortMetricName = ['Sorted by: ', SortMetric];
        end
        colormapStack_for_tiled(ID(i).graph.var_nfkb(ID(i).graph.order,:),ID(i).graph.celldata(ID(i).graph.order,:), ID(i).graph.opt_nfkb, fig_handle, axes.ax(i),SortMetricName); %calls subfunction that makes figure    
    end



    for i= 1:n
        axes.ax(i+n) = nexttile;
        if i < n
            ID(i).graph.opt_ktr.Name = '';
            SortMetricName = '';
        else
            SortMetricName = ['Sorted by: ', SortMetric];
        end
        colormapStack_for_tiled(ID(i).graph.var_ktr(ID(i).graph.order,:),ID(i).graph.celldata(ID(i).graph.order,:), ID(i).graph.opt_ktr, fig_handle,axes.ax(i+n), SortMetricName); %calls subfunction that makes figure    
    end
else
    figs.a = figure('name', 'HM NFkB');

    set(figs.a,'Position', [500 400 1000 300])
    fig_handle = gcf;
    tiledlayout(1,n)

    for i= 1:n
        axes.ax(i) = nexttile;
        if i < n
            ID(i).graph.opt_nfkb.Name = '';
            SortMetricName = '';
        else
            SortMetricName = ['Sorted by: ', SortMetric];
        end
        colormapStack_for_tiled(ID(i).graph.var_nfkb(ID(i).graph.order,:),ID(i).graph.celldata(ID(i).graph.order,:), ID(i).graph.opt_nfkb, fig_handle, axes.ax(i),SortMetricName); %calls subfunction that makes figure    
    end
end


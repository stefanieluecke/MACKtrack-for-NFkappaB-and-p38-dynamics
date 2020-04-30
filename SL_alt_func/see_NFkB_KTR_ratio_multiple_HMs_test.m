function [] = see_NFkB_KTR_ratio_multiple_HMs_test(id1,id2, id3, id4, id5, id6, varargin)

% 
%Enter highest responder first
%% INPUT PARSING
% Create input parser object, add required params from function input
p = inputParser;
% Required: ID input
valid_id = @(x) assert((isnumeric(x)&&length(x)==1)||exist(x,'file'),...
    'ID input must be spreadsheet ID or full file path');
addRequired(p,'id1',valid_id);

%Test 11/19/19 test using 6 ids as required parameters
%{
% Optional parameters
addParameter(p, 'id2', valid_id);
addParameter(p, 'id3', valid_id);
addParameter(p, 'id4', valid_id);
addParameter(p, 'id5', valid_id);
addParameter(p, 'id6', valid_id);
%}

% More required parameters
addRequired(p,'id2',valid_id);
addRequired(p,'id3',valid_id);
addRequired(p,'id4',valid_id);
addRequired(p,'id5',valid_id);
addRequired(p,'id6',valid_id);

expectedFlags = {'on','off'};
addParameter(p,'Verbose','on', @(x) any(validatestring(x,expectedFlags)));%checks whether optional name-value argument matches on or off %checks if x matches expectedFlags
valid_conv = @(x) assert(isnumeric(x)&&(x>=0)&&(length(x)==1),...
    'Parameter must be single integer >= 0'); %checks whether parameters below are single integers
addParameter(p,'ConvectionShift',1, valid_conv); %allows adjustment of convection shift (?)
addParameter(p,'MinLifetime',117, @isnumeric); %allows adjustment of minimum lifetime (?)
addParameter(p,'MinSize',90, valid_conv); %allows adjustment of minimum size (?)
addParameter(p,'TrimFrame',157, @isnumeric);
addParameter(p,'StartThreshNFkB',1000, valid_conv); %max allowable starting threshhold to filter out cells with pre-activated NFkB, default is 2
addParameter (p, 'OnThreshNFkB', 0, @isnumeric); %? not used in code?
addParameter (p, 'GraphLimitsNFkB',[-0.25 8],@isnumeric);
addParameter(p,'StartThreshKTR',1800, valid_conv); %max allowable starting threshhold to filter out cells with pre-activated KTR, default is 0.6
addParameter (p, 'OnThreshKTR', 0, @isnumeric);%? not used in code?
addParameter (p, 'GraphLimitsKTR',[-0.02 0.25],@isnumeric);
addParameter(p, 'SortMetric', 'peakfreq_nfkb');
expectedOrder = {'ascend', 'descend'};
addParameter(p, 'SortOrder', 'descend', @(x)any(validatestring(x, expectedOrder))); 
addParameter(p, 'StartTimePoint', 13, @isnumeric)

parse(p,id1,id2, id3, id4, id5, id6, varargin{:})
%%
[metrics1,aux1, graph1, info1, measure1] = nfkb_ktr_ratio_metrics_test(id1, 'MinLifetime',p.Results.MinLifetime,...
                            'ConvectionShift',p.Results.ConvectionShift, 'OnThreshNFkB',p.Results.OnThreshNFkB,'OnThreshKTR',p.Results.OnThreshKTR,...
                            'MinSize', p.Results.MinSize,'StartThreshNFkB', p.Results.StartThreshNFkB,'StartThreshKTR', p.Results.StartThreshKTR, 'Verbose', ... 
                            p.Results.Verbose, 'GraphLimitsNFkB', p.Results.GraphLimitsNFkB,'GraphLimitsKTR', p.Results.GraphLimitsKTR, 'TrimFrame', p.Results.TrimFrame, 'StartTimePoint', p.Results.StartTimePoint);

if isnumeric(p.Results.id2)
    [metrics2,aux2, graph2, info2, measure2] = nfkb_ktr_ratio_metrics_test(id2, 'MinLifetime',p.Results.MinLifetime,...
                            'ConvectionShift',p.Results.ConvectionShift, 'OnThreshNFkB',p.Results.OnThreshNFkB,'OnThreshKTR',p.Results.OnThreshKTR,...
                            'MinSize', p.Results.MinSize,'StartThreshNFkB', p.Results.StartThreshNFkB,'StartThreshKTR', p.Results.StartThreshKTR, 'Verbose', ... 
                            p.Results.Verbose, 'GraphLimitsNFkB', p.Results.GraphLimitsNFkB,'GraphLimitsKTR', p.Results.GraphLimitsKTR, 'TrimFrame', p.Results.TrimFrame, 'StartTimePoint', p.Results.StartTimePoint);
end 
if isnumeric(p.Results.id3)
    [metrics3,aux3, graph3, info3, measure3] = nfkb_ktr_ratio_metrics_test(id3, 'MinLifetime',p.Results.MinLifetime,...
                            'ConvectionShift',p.Results.ConvectionShift, 'OnThreshNFkB',p.Results.OnThreshNFkB,'OnThreshKTR',p.Results.OnThreshKTR,...
                            'MinSize', p.Results.MinSize,'StartThreshNFkB', p.Results.StartThreshNFkB,'StartThreshKTR', p.Results.StartThreshKTR, 'Verbose', ... 
                            p.Results.Verbose, 'GraphLimitsNFkB', p.Results.GraphLimitsNFkB,'GraphLimitsKTR', p.Results.GraphLimitsKTR, 'TrimFrame', p.Results.TrimFrame, 'StartTimePoint', p.Results.StartTimePoint);
end
if isnumeric(p.Results.id4)
    [metrics4,aux4, graph4, info4, measure4] = nfkb_ktr_ratio_metrics_test(id4, 'MinLifetime',p.Results.MinLifetime,...
                            'ConvectionShift',p.Results.ConvectionShift, 'OnThreshNFkB',p.Results.OnThreshNFkB,'OnThreshKTR',p.Results.OnThreshKTR,...
                            'MinSize', p.Results.MinSize,'StartThreshNFkB', p.Results.StartThreshNFkB,'StartThreshKTR', p.Results.StartThreshKTR, 'Verbose', ... 
                            p.Results.Verbose, 'GraphLimitsNFkB', p.Results.GraphLimitsNFkB,'GraphLimitsKTR', p.Results.GraphLimitsKTR, 'TrimFrame', p.Results.TrimFrame, 'StartTimePoint', p.Results.StartTimePoint);
end
if isnumeric(p.Results.id5)
    [metrics5,aux5, graph5, info5, measure5] = nfkb_ktr_ratio_metrics_test(id5, 'MinLifetime',p.Results.MinLifetime,...
                            'ConvectionShift',p.Results.ConvectionShift, 'OnThreshNFkB',p.Results.OnThreshNFkB,'OnThreshKTR',p.Results.OnThreshKTR,...
                            'MinSize', p.Results.MinSize,'StartThreshNFkB', p.Results.StartThreshNFkB,'StartThreshKTR', p.Results.StartThreshKTR, 'Verbose', ... 
                            p.Results.Verbose, 'GraphLimitsNFkB', p.Results.GraphLimitsNFkB,'GraphLimitsKTR', p.Results.GraphLimitsKTR, 'TrimFrame', p.Results.TrimFrame, 'StartTimePoint', p.Results.StartTimePoint);
end
if isnumeric(p.Results.id6)
    [metrics6,aux6, graph6, info6, measure6] = nfkb_ktr_ratio_metrics_test(id6, 'MinLifetime',p.Results.MinLifetime,...
                            'ConvectionShift',p.Results.ConvectionShift, 'OnThreshNFkB',p.Results.OnThreshNFkB,'OnThreshKTR',p.Results.OnThreshKTR,...
                            'MinSize', p.Results.MinSize,'StartThreshNFkB', p.Results.StartThreshNFkB,'StartThreshKTR', p.Results.StartThreshKTR, 'Verbose', ... 
                            p.Results.Verbose, 'GraphLimitsNFkB', p.Results.GraphLimitsNFkB,'GraphLimitsKTR', p.Results.GraphLimitsKTR, 'TrimFrame', p.Results.TrimFrame, 'StartTimePoint', p.Results.StartTimePoint);
end


%%

%{
SortMetric = p.Results.SortMetric;

[~,graph1.order] = sort(metrics1.(SortMetric), p.Results.SortOrder);
figs.a = figure('name', 'Heatmaps');
set(figs.a,'Position', [500 400 1000 600])
colormapStack_double(graph1.var_nfkb(graph1.order,:), graph1.var_ktr(graph1.order,:),graph1.celldata(graph1.order,:), graph1.opt_nfkb, graph1.opt_ktr);


if isnumeric(p.Results.id2) 
    [~,graph2.order] = sort(metrics2.(SortMetric), p.Results.SortOrder);
    figs.b = figure('name', 'Heatmaps');
    set(figs.b,'Position', [500 400 1000 600])
    colormapStack_double(graph2.var_nfkb(graph2.order,:), graph2.var_ktr(graph2.order,:),graph2.celldata(graph2.order,:), graph2.opt_nfkb, graph2.opt_ktr);

    
    figs.c = figure;
    sub_a = subplot(1,2,1);
    sub_b = subplot(1,2,2);
   
    copyobj(allchild(get(figs.a,'CurrentAxes')),sub_a);
    copyobj(allchild(get(figs.b,'CurrentAxes')),sub_b);
    
    
    copyobj(figs.a.Children,sub_a.Children);
    copyobj(figs.b.Children,sub_b.Children);
    
    %get_a = get(figs.a);
    %get_b = get(figs.b);
    %set(sub_a,get_a);
    %set(sub_b,get_b);
    %copyobj(get(figs.a),sub_a);
    %copyobj(get(figs.b),sub_b);

end    
 %}

%% Attempt to make NFkB and KTR HMs separately
%{
SortMetric = p.Results.SortMetric;

[~,graph1.order] = sort(metrics1.(SortMetric), p.Results.SortOrder);

figs.a = figure('name','Heatmap_NFkB');
set(figs.a,'Position', [500 7 400 600])
subplot(1,6,1);
colormapStack(graph1.var_nfkb(graph1.order,:),graph1.celldata(graph1.order,:), graph1.opt_nfkb); %calls subfunction that makes figure    

if isnumeric(p.Results.id2) 
    [~,graph2.order] = sort(metrics2.(SortMetric), p.Results.SortOrder);
    subplot(1,6,2);
    %figs.b = figure('name', 'Heatmaps');
    %set(figs.b,'Position', [500 400 1000 600])
    colormapStack(graph2.var_nfkb(graph2.order,:),graph2.celldata(graph2.order,:), graph2.opt_nfkb);
end
%} 
%% Attempt to hardcode this in colormapStack_multiple

SortMetric = p.Results.SortMetric;
[~,graph1.order] = sort(metrics1.(SortMetric), p.Results.SortOrder);
[~,graph2.order] = sort(metrics2.(SortMetric), p.Results.SortOrder);
[~,graph3.order] = sort(metrics3.(SortMetric), p.Results.SortOrder);
[~,graph4.order] = sort(metrics4.(SortMetric), p.Results.SortOrder);
[~,graph5.order] = sort(metrics5.(SortMetric), p.Results.SortOrder);
[~,graph6.order] = sort(metrics6.(SortMetric), p.Results.SortOrder);

% Heatmap
HMinput1 = struct('var_nfkb_sorted',graph1.var_nfkb(graph1.order,:), 'var_ktr_sorted', graph1.var_ktr(graph1.order,:),...
    'celldata',graph1.celldata(graph1.order,:), 'opt_nfkb', graph1.opt_nfkb, 'opt_ktr', graph1.opt_ktr, 'title', info1.name);
HMinput2 = struct('var_nfkb_sorted',graph2.var_nfkb(graph2.order,:), 'var_ktr_sorted', graph2.var_ktr(graph2.order,:),...
    'celldata',graph2.celldata(graph2.order,:), 'opt_nfkb', graph2.opt_nfkb, 'opt_ktr', graph2.opt_ktr, 'title', info2.name);
HMinput3 = struct('var_nfkb_sorted',graph3.var_nfkb(graph3.order,:), 'var_ktr_sorted', graph3.var_ktr(graph3.order,:),...
    'celldata',graph3.celldata(graph3.order,:), 'opt_nfkb', graph3.opt_nfkb, 'opt_ktr', graph3.opt_ktr, 'title', info3.name);
HMinput4 = struct('var_nfkb_sorted',graph4.var_nfkb(graph4.order,:), 'var_ktr_sorted', graph4.var_ktr(graph4.order,:),...
    'celldata',graph4.celldata(graph4.order,:), 'opt_nfkb', graph4.opt_nfkb, 'opt_ktr', graph4.opt_ktr, 'title', info4.name);
HMinput5 = struct('var_nfkb_sorted',graph5.var_nfkb(graph5.order,:), 'var_ktr_sorted', graph5.var_ktr(graph5.order,:),...
    'celldata',graph5.celldata(graph5.order,:), 'opt_nfkb', graph5.opt_nfkb, 'opt_ktr', graph5.opt_ktr, 'title', info5.name);
HMinput6 = struct('var_nfkb_sorted',graph6.var_nfkb(graph6.order,:), 'var_ktr_sorted', graph6.var_ktr(graph6.order,:),...
    'celldata',graph6.celldata(graph6.order,:), 'opt_nfkb', graph6.opt_nfkb, 'opt_ktr', graph6.opt_ktr, 'title', info6.name);


% Heatmap NFkB and KTR plots next to each other sorted the same
figs.e = figure('name', 'Multiple_Heatmaps_NFkB+KTR');
set(figs.e,'Position', [50 100 2000 900])
colormapStack_multiple_double(HMinput1,HMinput2,HMinput3,HMinput4,HMinput5,HMinput6);

%{
f=6;
for i = 1:f
    %fname = [num2str(i)];
    %HMinput(fname) = struct('var_nfkb_sorted', graph(fname).var_nfkb(graph(fname).order,:), 'var_ktr_sorted', graph(fname).var_ktr(graph(fname).order,:),graph(fname).celldata(graph(fname).order,:), graph(fname).opt_nfkb, graph(fname).opt_ktr);
    HMinput(i) = struct('var_nfkb_sorted', grap(i).var_nfkb(grap(i).order,:), 'var_ktr_sorted', grap(i).var_ktr(grap(i).order,:),....
        'celldata', grap(i).celldata(grap(i).order,:), 'opt_nfkb', grap(i).opt_nfkb, 'opt_ktr', grap(i).opt_ktr);
end
%}

%Heatmaps NFkB plots next to each other
figs.f = figure('name', 'Multiple_Heatmaps_NFkB');
set(figs.f,'Position', [50 100 1800 500])
colormapStack_multiple_NFkB(HMinput1,HMinput2,HMinput3,HMinput4,HMinput5,HMinput6);


%Heatmaps KTR plots next to each other
figs.f = figure('name', 'Multiple_Heatmaps_KTR');
set(figs.f,'Position', [50 100 1800 500])
colormapStack_multiple_KTR(HMinput1,HMinput2,HMinput3,HMinput4,HMinput5,HMinput6);

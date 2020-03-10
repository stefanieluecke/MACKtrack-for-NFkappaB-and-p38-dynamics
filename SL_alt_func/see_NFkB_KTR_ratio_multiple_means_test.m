function [] = see_NFkB_KTR_ratio_multiple_means_test(id1, varargin)

% 
%Enter highest responder first
%% INPUT PARSING
% Create input parser object, add required params from function input
p = inputParser;
% Required: ID input
valid_id = @(x) assert((isnumeric(x)&&length(x)==1)||exist(x,'file'),...
    'ID input must be spreadsheet ID or full file path');
addRequired(p,'id1',valid_id);

% Optional parameters
addParameter(p, 'id2', valid_id);
addParameter(p, 'id3', valid_id);
addParameter(p, 'id4', valid_id);
addParameter(p, 'id5', valid_id);
addParameter(p, 'id6', valid_id);

expectedFlags = {'on','off'};
addParameter(p,'Verbose','on', @(x) any(validatestring(x,expectedFlags)));%checks whether optional name-value argument matches on or off %checks if x matches expectedFlags
valid_conv = @(x) assert(isnumeric(x)&&(x>=0)&&(length(x)==1),...
    'Parameter must be single integer >= 0'); %checks whether parameters below are single integers
addParameter(p,'ConvectionShift',1, valid_conv); %allows adjustment of convection shift (?)
addParameter(p,'MinLifetime',117, @isnumeric); %allows adjustment of minimum lifetime (?)
addParameter(p,'MinSize',90, valid_conv); %allows adjustment of minimum size (?)
addParameter(p,'TrimFrame',157, @isnumeric);
addParameter(p,'StartThreshNFkB',1, valid_conv); %max allowable starting threshhold to filter out cells with pre-activated NFkB, default is 2
addParameter (p, 'OnThreshNFkB', 0, @isnumeric); %? not used in code?
addParameter (p, 'GraphLimitsNFkB',[-0.25 8],@isnumeric);
addParameter(p,'StartThreshKTR',1, valid_conv); %max allowable starting threshhold to filter out cells with pre-activated KTR, default is 0.6
addParameter (p, 'OnThreshKTR', 0, @isnumeric);%? not used in code?
addParameter (p, 'GraphLimitsKTR',[0 500],@isnumeric);
addParameter(p, 'SortMetric', 'peakfreq_nfkb');
expectedOrder = {'ascend', 'descend'};
addParameter(p, 'SortOrder', 'descend', @(x)any(validatestring(x, expectedOrder))); 
addParameter(p, 'StartTimePoint', 13, @isnumeric)

parse(p,id1, varargin{:})
%%
[metrics1,aux1, graph1, info1, measure1] = nfkb_ktr_ratio_metrics_test(id1, 'MinLifetime',p.Results.MinLifetime,...
                            'ConvectionShift',p.Results.ConvectionShift, 'OnThreshNFkB',p.Results.OnThreshNFkB,'OnThreshKTR',p.Results.OnThreshKTR,...
                            'MinSize', p.Results.MinSize,'StartThreshNFkB', p.Results.StartThreshNFkB,'StartThreshKTR', p.Results.StartThreshKTR, 'Verbose', ... 
                            p.Results.Verbose, 'GraphLimitsNFkB', p.Results.GraphLimitsNFkB,'GraphLimitsKTR', p.Results.GraphLimitsKTR, 'TrimFrame', p.Results.TrimFrame, 'StartTimePoint', p.Results.StartTimePoint);

if isnumeric(p.Results.id2)
    [metrics2,aux2, graph2, info2, measure2] = nfkb_ktr_ratio_metrics_test(p.Results.id2, 'MinLifetime',p.Results.MinLifetime,...
                            'ConvectionShift',p.Results.ConvectionShift, 'OnThreshNFkB',p.Results.OnThreshNFkB,'OnThreshKTR',p.Results.OnThreshKTR,...
                            'MinSize', p.Results.MinSize,'StartThreshNFkB', p.Results.StartThreshNFkB,'StartThreshKTR', p.Results.StartThreshKTR, 'Verbose', ... 
                            p.Results.Verbose, 'GraphLimitsNFkB', p.Results.GraphLimitsNFkB,'GraphLimitsKTR', p.Results.GraphLimitsKTR, 'TrimFrame', p.Results.TrimFrame, 'StartTimePoint', p.Results.StartTimePoint);
end 
if isnumeric(p.Results.id3)
    [metrics3,aux3, graph3, info3, measure3] = nfkb_ktr_ratio_metrics_test(p.Results.id3, 'MinLifetime',p.Results.MinLifetime,...
                            'ConvectionShift',p.Results.ConvectionShift, 'OnThreshNFkB',p.Results.OnThreshNFkB,'OnThreshKTR',p.Results.OnThreshKTR,...
                            'MinSize', p.Results.MinSize,'StartThreshNFkB', p.Results.StartThreshNFkB,'StartThreshKTR', p.Results.StartThreshKTR, 'Verbose', ... 
                            p.Results.Verbose, 'GraphLimitsNFkB', p.Results.GraphLimitsNFkB,'GraphLimitsKTR', p.Results.GraphLimitsKTR, 'TrimFrame', p.Results.TrimFrame, 'StartTimePoint', p.Results.StartTimePoint);
end
if isnumeric(p.Results.id4)
    [metrics4,aux4, graph4, info4, measure4] = nfkb_ktr_ratio_metrics_test(p.Results.id4, 'MinLifetime',p.Results.MinLifetime,...
                            'ConvectionShift',p.Results.ConvectionShift, 'OnThreshNFkB',p.Results.OnThreshNFkB,'OnThreshKTR',p.Results.OnThreshKTR,...
                            'MinSize', p.Results.MinSize,'StartThreshNFkB', p.Results.StartThreshNFkB,'StartThreshKTR', p.Results.StartThreshKTR, 'Verbose', ... 
                            p.Results.Verbose, 'GraphLimitsNFkB', p.Results.GraphLimitsNFkB,'GraphLimitsKTR', p.Results.GraphLimitsKTR, 'TrimFrame', p.Results.TrimFrame, 'StartTimePoint', p.Results.StartTimePoint);
end
if isnumeric(p.Results.id5)
    [metrics5,aux5, graph5, info5, measure5] = nfkb_ktr_ratio_metrics_test(p.Results.id5, 'MinLifetime',p.Results.MinLifetime,...
                            'ConvectionShift',p.Results.ConvectionShift, 'OnThreshNFkB',p.Results.OnThreshNFkB,'OnThreshKTR',p.Results.OnThreshKTR,...
                            'MinSize', p.Results.MinSize,'StartThreshNFkB', p.Results.StartThreshNFkB,'StartThreshKTR', p.Results.StartThreshKTR, 'Verbose', ... 
                            p.Results.Verbose, 'GraphLimitsNFkB', p.Results.GraphLimitsNFkB,'GraphLimitsKTR', p.Results.GraphLimitsKTR, 'TrimFrame', p.Results.TrimFrame, 'StartTimePoint', p.Results.StartTimePoint);
end
if isnumeric(p.Results.id6)
    [metrics6,aux6, graph6, info6, measure6] = nfkb_ktr_ratio_metrics_test(p.Results.id6, 'MinLifetime',p.Results.MinLifetime,...
                            'ConvectionShift',p.Results.ConvectionShift, 'OnThreshNFkB',p.Results.OnThreshNFkB,'OnThreshKTR',p.Results.OnThreshKTR,...
                            'MinSize', p.Results.MinSize,'StartThreshNFkB', p.Results.StartThreshNFkB,'StartThreshKTR', p.Results.StartThreshKTR, 'Verbose', ... 
                            p.Results.Verbose, 'GraphLimitsNFkB', p.Results.GraphLimitsNFkB,'GraphLimitsKTR', p.Results.GraphLimitsKTR, 'TrimFrame', p.Results.TrimFrame, 'StartTimePoint', p.Results.StartTimePoint);
end

%{
for i=2:6
    [metrics(i),aux, graph, info, measure(i)] = nfkb_ktr_ratio_metrics(id1, 'MinLifetime',p.Results.MinLifetime,...
                            'ConvectionShift',p.Results.ConvectionShift, 'OnThreshNFkB',p.Results.OnThreshNFkB,'OnThreshKTR',p.Results.OnThreshKTR,...
                            'MinSize', p.Results.MinSize,'StartThreshNFkB', p.Results.StartThreshNFkB,'StartThreshKTR', p.Results.StartThreshKTR, 'Verbose', ... 
                            p.Results.Verbose, 'GraphLimitsNFkB', p.Results.GraphLimitsNFkB,'GraphLimitsKTR', p.Results.GraphLimitsKTR, 'TrimFrame', p.Results.TrimFrame, 'StartTimePoint', p.Results.StartTimePoint);
%}                        
%%
%add in line plot of averages 
    colors = setcolors;    
% Line plot NFkB(mean+/-std)
    graph1.line.top_nfkb = nanmean(graph1.var_nfkb) + nanstd(graph1.var_nfkb); %creates upper bound
    graph1.line.bot_nfkb = nanmean(graph1.var_nfkb) - nanstd(graph1.var_nfkb); %creates lower bound
     
    figs.c = figure('name','OverlayMeanTrajectoriesNFkB'); %creates figure and names it
    set(figs.c,'Position', [555   743   1145   500]) %positions figure within window
    c = fill([graph1.t,graph1.t(end:-1:1)],[graph1.line.top_nfkb,graph1.line.bot_nfkb(end:-1:1)],colors.dark_red); %creates filled polygon determined by time axis and upper and lower bound
    set(c, 'EdgeColor', colors.dark_red);
    alpha(c,0.3)%adds transparency
    hold on
    graph1.line.main_nfkb = plot(graph1.t, nanmean(graph1.var_nfkb)); %adds the mean to the figure
    set(graph1.line.main_nfkb,'Color', colors.dark_red,'LineWidth',3); %determined properties of the mean line
    hold on
    legend([graph1.line.main_nfkb],{info1.name}, 'Interpreter', 'none')
    
    set(gca,'XTick',graph1.opt_nfkb.TimeTicks,'YTick',graph1.opt_nfkb.MeasurementTicks) %sets axis tick properites
    set(gca,'YTickLabel',graph1.opt_nfkb.MeasurementTickLabels,'TickLength',[0.005 0.005]) %sets tick label properties
   
    ylabel('NFkB activation', 'FontSize', 14);
    xlabel('Time (h)','FontSize',14);
    axis([min(graph1.opt_nfkb.Times) max(graph1.opt_nfkb.Times) graph1.opt_nfkb.MeasurementBounds]) %sets axis to determined timescale and determined measurement bounds, derived from maketicks function

if isnumeric(p.Results.id2)    
    graph2.line.top_nfkb = nanmean(graph2.var_nfkb) + nanstd(graph2.var_nfkb); %creates upper bound
    graph2.line.bot_nfkb = nanmean(graph2.var_nfkb) - nanstd(graph2.var_nfkb); %creates lower bound
    c = fill([graph2.t,graph2.t(end:-1:1)],[graph2.line.top_nfkb,graph2.line.bot_nfkb(end:-1:1)],colors.light_red); %creates filled polygon determined by time axis and upper and lower bound
    set(c, 'EdgeColor', colors.light_red);
    alpha(c,0.3)%adds transparency
    hold on
    graph2.line.main_nfkb = plot(graph2.t, nanmean(graph2.var_nfkb)); %adds the mean to the figure
    set(graph2.line.main_nfkb,'Color',colors.light_red,'LineWidth',3); %determined properties of the mean line
      legend([graph1.line.main_nfkb,graph2.line.main_nfkb],{info1.name,info2.name}, 'Interpreter', 'none')
    hold on
     
end

if isnumeric(p.Results.id3)    
    graph3.line.top_nfkb = nanmean(graph3.var_nfkb) + nanstd(graph3.var_nfkb); %creates upper bound
    graph3.line.bot_nfkb = nanmean(graph3.var_nfkb) - nanstd(graph3.var_nfkb); %creates lower bound
    c = fill([graph3.t,graph3.t(end:-1:1)],[graph3.line.top_nfkb,graph3.line.bot_nfkb(end:-1:1)],colors.orange2); %creates filled polygon determined by time axis and upper and lower bound
    set(c, 'EdgeColor', colors.orange2);
    alpha(c,0.3)%adds transparency
    hold on
    graph3.line.main_nfkb = plot(graph3.t, nanmean(graph3.var_nfkb)); %adds the mean to the figure
    set(graph3.line.main_nfkb,'Color',colors.orange2,'LineWidth',3); %determined properties of the mean line
    hold on
    legend([graph1.line.main_nfkb,graph2.line.main_nfkb,graph3.line.main_nfkb],{info1.name,info2.name,info3.name}, 'Interpreter', 'none')

end
    if isnumeric(p.Results.id4)    
    graph4.line.top_nfkb = nanmean(graph4.var_nfkb) + nanstd(graph4.var_nfkb); %creates upper bound
    graph4.line.bot_nfkb = nanmean(graph4.var_nfkb) - nanstd(graph4.var_nfkb); %creates lower bound
    c = fill([graph4.t,graph4.t(end:-1:1)],[graph4.line.top_nfkb,graph4.line.bot_nfkb(end:-1:1)],colors.turquoise); %creates filled polygon determined by time axis and upper and lower bound
    set(c, 'EdgeColor', colors.turquoise);
    alpha(c,0.3)%adds transparency
    hold on
    graph4.line.main_nfkb = plot(graph4.t, nanmean(graph4.var_nfkb)); %adds the mean to the figure
    set(graph4.line.main_nfkb,'Color',colors.turquoise,'LineWidth',3); %determined properties of the mean line
    hold on
    legend([graph1.line.main_nfkb,graph2.line.main_nfkb,graph3.line.main_nfkb,graph4.line.main_nfkb],{info1.name,info2.name,info3.name,info4.name}, 'Interpreter', 'none')

end
if isnumeric(p.Results.id5)    
    graph5.line.top_nfkb = nanmean(graph5.var_nfkb) + nanstd(graph5.var_nfkb); %creates upper bound
    graph5.line.bot_nfkb = nanmean(graph5.var_nfkb) - nanstd(graph5.var_nfkb); %creates lower bound
    c = fill([graph5.t,graph5.t(end:-1:1)],[graph5.line.top_nfkb,graph5.line.bot_nfkb(end:-1:1)],colors.sky_blue); %creates filled polygon determined by time axis and upper and lower bound
        set(c, 'EdgeColor', colors.sky_blue);
    alpha(c,0.3)%adds transparency
    hold on
    graph5.line.main_nfkb = plot(graph5.t, nanmean(graph5.var_nfkb)); %adds the mean to the figure
    set(graph5.line.main_nfkb,'Color',colors.sky_blue,'LineWidth',3); %determined properties of the mean line
    hold on
    legend([graph1.line.main_nfkb,graph2.line.main_nfkb,graph3.line.main_nfkb,graph4.line.main_nfkb,graph5.line.main_nfkb],{info1.name,info2.name,info3.name,info4.name,info5.name}, 'Interpreter', 'none')
    
end
if isnumeric(p.Results.id6)    
    graph6.line.top_nfkb = nanmean(graph6.var_nfkb) + nanstd(graph6.var_nfkb); %creates upper bound
    graph6.line.bot_nfkb = nanmean(graph6.var_nfkb) - nanstd(graph6.var_nfkb); %creates lower bound
    c = fill([graph6.t,graph6.t(end:-1:1)],[graph6.line.top_nfkb,graph6.line.bot_nfkb(end:-1:1)],colors.blue); %creates filled polygon determined by time axis and upper and lower bound
        set(c, 'EdgeColor', colors.blue);
    alpha(c,0.3)%adds transparency
    hold on
    graph6.line.main_nfkb = plot(graph6.t, nanmean(graph6.var_nfkb)); %adds the mean to the figure
    set(graph6.line.main_nfkb,'Color',colors.blue,'LineWidth',3); %determined properties of the mean line
    % legend(info1.name,info2.name,info3.name,info4.name,info5.name,info6.name, 'Interpreter', 'none');
    legend([graph1.line.main_nfkb,graph2.line.main_nfkb,graph3.line.main_nfkb,graph4.line.main_nfkb,graph5.line.main_nfkb,graph6.line.main_nfkb],{info1.name,info2.name,info3.name,info4.name,info5.name,info6.name}, 'Interpreter', 'none')
    hold off
    
end
    
%%
% Line plot KTR(mean+/-std)
    graph1.line.top_ktr = nanmean(graph1.var_ktr) + nanstd(graph1.var_ktr); %creates upper bound
    graph1.line.bot_ktr = nanmean(graph1.var_ktr) - nanstd(graph1.var_ktr); %creates lower bound
    
    figs.d = figure('name','OverlayMeanTrajectoriesKTR'); %creates figure and names it
    set(figs.d,'Position', [555   743   1145   500]) %positions figure within window
    c = fill([graph1.t,graph1.t(end:-1:1)],[graph1.line.top_ktr,graph1.line.bot_ktr(end:-1:1)],colors.dark_red); %creates filled polygon determined by time axis and upper and lower bound
        set(c, 'EdgeColor', colors.dark_red);
    alpha(c,0.4)%adds transparency
    hold on
    graph1.line.main_ktr = plot(graph1.t, nanmean(graph1.var_ktr)); %adds the mean to the figure
    set(graph1.line.main_ktr,'Color', colors.dark_red,'LineWidth',3); %determined properties of the mean line
    hold on
     legend([graph1.line.main_ktr],{info1.name}, 'Interpreter', 'none');
    
    set(gca,'XTick',graph1.opt_ktr.TimeTicks,'YTick',graph1.opt_ktr.MeasurementTicks) %sets axis tick properites
    set(gca,'YTickLabel',graph1.opt_ktr.MeasurementTickLabels,'TickLength',[0.005 0.005]) %sets tick label properties
    ylabel('KTR activation', 'FontSize', 14);
    xlabel('Time (h)','FontSize',14);
    axis([min(graph1.opt_ktr.Times) max(graph1.opt_ktr.Times) graph1.opt_ktr.MeasurementBounds]) %sets axis to determined timescale and determined measurement bounds, derived from maketicks function

if isnumeric(p.Results.id2)    
    graph2.line.top_ktr = nanmean(graph2.var_ktr) + nanstd(graph2.var_ktr); %creates upper bound
    graph2.line.bot_ktr = nanmean(graph2.var_ktr) - nanstd(graph2.var_ktr); %creates lower bound
    c = fill([graph2.t,graph2.t(end:-1:1)],[graph2.line.top_ktr,graph2.line.bot_ktr(end:-1:1)],colors.light_red); %creates filled polygon determined by time axis and upper and lower bound
    set(c, 'EdgeColor', colors.light_red);
    alpha(c,0.5)%adds transparency
    hold on
    graph2.line.main_ktr = plot(graph2.t, nanmean(graph2.var_ktr)); %adds the mean to the figure
    set(graph2.line.main_ktr,'Color',colors.light_red,'LineWidth',3); %determined properties of the mean line
    hold on
    legend([graph1.line.main_ktr,graph2.line.main_ktr],{info1.name,info2.name}, 'Interpreter', 'none')
end

if isnumeric(p.Results.id3)    
    graph3.line.top_ktr = nanmean(graph3.var_ktr) + nanstd(graph3.var_ktr); %creates upper bound
    graph3.line.bot_ktr = nanmean(graph3.var_ktr) - nanstd(graph3.var_ktr); %creates lower bound
    c = fill([graph3.t,graph3.t(end:-1:1)],[graph3.line.top_ktr,graph3.line.bot_ktr(end:-1:1)],colors.orange2); %creates filled polygon determined by time axis and upper and lower bound
    set(c, 'EdgeColor', colors.orange2);
    alpha(c,0.3)%adds transparency
    hold on
    graph3.line.main_ktr = plot(graph3.t, nanmean(graph3.var_ktr)); %adds the mean to the figure
    set(graph3.line.main_ktr,'Color',colors.orange2,'LineWidth',3); %determined properties of the mean line
    hold on
    legend([graph1.line.main_ktr,graph2.line.main_ktr,graph3.line.main_ktr],{info1.name,info2.name,info3.name}, 'Interpreter', 'none')
end
    if isnumeric(p.Results.id4)    
    graph4.line.top_ktr = nanmean(graph4.var_ktr) + nanstd(graph4.var_ktr); %creates upper bound
    graph4.line.bot_ktr = nanmean(graph4.var_ktr) - nanstd(graph4.var_ktr); %creates lower bound
    c = fill([graph4.t,graph4.t(end:-1:1)],[graph4.line.top_ktr,graph4.line.bot_ktr(end:-1:1)],colors.turquoise); %creates filled polygon determined by time axis and upper and lower bound
        set(c, 'EdgeColor', colors.turquoise);
    alpha(c,0.3)%adds transparency
    hold on
    graph4.line.main_ktr = plot(graph4.t, nanmean(graph4.var_ktr)); %adds the mean to the figure
    set(graph4.line.main_ktr,'Color',colors.turquoise,'LineWidth',3); %determined properties of the mean line
    hold on
legend([graph1.line.main_ktr,graph2.line.main_ktr,graph3.line.main_ktr,graph4.line.main_ktr],{info1.name,info2.name,info3.name,info4.name}, 'Interpreter', 'none')
end
if isnumeric(p.Results.id5)    
    graph5.line.top_ktr = nanmean(graph5.var_ktr) + nanstd(graph5.var_ktr); %creates upper bound
    graph5.line.bot_ktr = nanmean(graph5.var_ktr) - nanstd(graph5.var_ktr); %creates lower bound
    c = fill([graph5.t,graph5.t(end:-1:1)],[graph5.line.top_ktr,graph5.line.bot_ktr(end:-1:1)],colors.sky_blue); %creates filled polygon determined by time axis and upper and lower bound
        set(c, 'EdgeColor', colors.sky_blue);
    alpha(c,0.2)%adds transparency
    hold on
    graph5.line.main_ktr = plot(graph5.t, nanmean(graph5.var_ktr)); %adds the mean to the figure
    set(graph5.line.main_ktr,'Color',colors.sky_blue,'LineWidth',3); %determined properties of the mean line
    hold on
    legend([graph1.line.main_ktr,graph2.line.main_ktr,graph3.line.main_ktr,graph4.line.main_ktr,graph5.line.main_ktr],{info1.name,info2.name,info3.name,info4.name,info5.name}, 'Interpreter', 'none')
end
if isnumeric(p.Results.id6)    
    graph6.line.top_ktr = nanmean(graph6.var_ktr) + nanstd(graph6.var_ktr); %creates upper bound
    graph6.line.bot_ktr = nanmean(graph6.var_ktr) - nanstd(graph6.var_ktr); %creates lower bound
    c = fill([graph6.t,graph6.t(end:-1:1)],[graph6.line.top_ktr,graph6.line.bot_ktr(end:-1:1)],colors.blue); %creates filled polygon determined by time axis and upper and lower bound
        set(c, 'EdgeColor', colors.blue);
    alpha(c, 0.2)%adds transparency
    hold on
    graph6.line.main_ktr = plot(graph6.t, nanmean(graph6.var_ktr)); %adds the mean to the figure
    set(graph6.line.main_ktr,'Color',colors.blue,'LineWidth',3); %determined properties of the mean line
    legend([graph1.line.main_ktr,graph2.line.main_ktr,graph3.line.main_ktr,graph4.line.main_ktr,graph5.line.main_ktr,graph6.line.main_ktr],{info1.name,info2.name,info3.name,info4.name,info5.name,info6.name}, 'Interpreter', 'none')
 
    hold off
end

%%
% Line plot KTR and NFkB(mean+/-std)
    figs.e =figure('name', 'OverlayMeanTrajectories');
    set(figs.e,'Position', [555   743   1100   1000]) %positions figure within window
% NFkB plot
    subplot(2,1,1)   
    
    graph1.line.top_nfkb = nanmean(graph1.var_nfkb) + nanstd(graph1.var_nfkb); %creates upper bound
    graph1.line.bot_nfkb = nanmean(graph1.var_nfkb) - nanstd(graph1.var_nfkb); %creates lower bound
    c = fill([graph1.t,graph1.t(end:-1:1)],[graph1.line.top_nfkb,graph1.line.bot_nfkb(end:-1:1)],colors.dark_red); %creates filled polygon determined by time axis and upper and lower bound
    set(c, 'EdgeColor', colors.dark_red);
    alpha(c,0.3)%adds transparency
    hold on
    graph1.line.main_nfkb = plot(graph1.t, nanmean(graph1.var_nfkb)); %adds the mean to the figure
    set(graph1.line.main_nfkb,'Color', colors.dark_red,'LineWidth',3); %determined properties of the mean line
    hold on
    legend([graph1.line.main_nfkb],{info1.name}, 'Interpreter', 'none')
    
    set(gca,'XTick',graph1.opt_nfkb.TimeTicks,'YTick',graph1.opt_nfkb.MeasurementTicks) %sets axis tick properites
    set(gca,'YTickLabel',graph1.opt_nfkb.MeasurementTickLabels,'TickLength',[0.005 0.005]) %sets tick label properties
   
    ylabel('NFkB activation', 'FontSize', 14);
    xlabel('Time (h)','FontSize',14);
    axis([min(graph1.opt_nfkb.Times) max(graph1.opt_nfkb.Times) graph1.opt_nfkb.MeasurementBounds]) %sets axis to determined timescale and determined measurement bounds, derived from maketicks function

if isnumeric(p.Results.id2)    
    graph2.line.top_nfkb = nanmean(graph2.var_nfkb) + nanstd(graph2.var_nfkb); %creates upper bound
    graph2.line.bot_nfkb = nanmean(graph2.var_nfkb) - nanstd(graph2.var_nfkb); %creates lower bound
    c = fill([graph2.t,graph2.t(end:-1:1)],[graph2.line.top_nfkb,graph2.line.bot_nfkb(end:-1:1)],colors.light_red); %creates filled polygon determined by time axis and upper and lower bound
    set(c, 'EdgeColor', colors.light_red);
    alpha(c,0.3)%adds transparency
    hold on
    graph2.line.main_nfkb = plot(graph2.t, nanmean(graph2.var_nfkb)); %adds the mean to the figure
    set(graph2.line.main_nfkb,'Color',colors.light_red,'LineWidth',3); %determined properties of the mean line
      legend([graph1.line.main_nfkb,graph2.line.main_nfkb],{info1.name,info2.name}, 'Interpreter', 'none')
    hold on
     
end

if isnumeric(p.Results.id3)    
    graph3.line.top_nfkb = nanmean(graph3.var_nfkb) + nanstd(graph3.var_nfkb); %creates upper bound
    graph3.line.bot_nfkb = nanmean(graph3.var_nfkb) - nanstd(graph3.var_nfkb); %creates lower bound
    c = fill([graph3.t,graph3.t(end:-1:1)],[graph3.line.top_nfkb,graph3.line.bot_nfkb(end:-1:1)],colors.orange2); %creates filled polygon determined by time axis and upper and lower bound
    set(c, 'EdgeColor', colors.orange2);
    alpha(c,0.3)%adds transparency
    hold on
    graph3.line.main_nfkb = plot(graph3.t, nanmean(graph3.var_nfkb)); %adds the mean to the figure
    set(graph3.line.main_nfkb,'Color',colors.orange2,'LineWidth',3); %determined properties of the mean line
    hold on
    legend([graph1.line.main_nfkb,graph2.line.main_nfkb,graph3.line.main_nfkb],{info1.name,info2.name,info3.name}, 'Interpreter', 'none')

end
    if isnumeric(p.Results.id4)    
    graph4.line.top_nfkb = nanmean(graph4.var_nfkb) + nanstd(graph4.var_nfkb); %creates upper bound
    graph4.line.bot_nfkb = nanmean(graph4.var_nfkb) - nanstd(graph4.var_nfkb); %creates lower bound
    c = fill([graph4.t,graph4.t(end:-1:1)],[graph4.line.top_nfkb,graph4.line.bot_nfkb(end:-1:1)],colors.turquoise); %creates filled polygon determined by time axis and upper and lower bound
    set(c, 'EdgeColor', colors.turquoise);
    alpha(c,0.3)%adds transparency
    hold on
    graph4.line.main_nfkb = plot(graph4.t, nanmean(graph4.var_nfkb)); %adds the mean to the figure
    set(graph4.line.main_nfkb,'Color',colors.turquoise,'LineWidth',3); %determined properties of the mean line
    hold on
    legend([graph1.line.main_nfkb,graph2.line.main_nfkb,graph3.line.main_nfkb,graph4.line.main_nfkb],{info1.name,info2.name,info3.name,info4.name}, 'Interpreter', 'none')

end
if isnumeric(p.Results.id5)    
    graph5.line.top_nfkb = nanmean(graph5.var_nfkb) + nanstd(graph5.var_nfkb); %creates upper bound
    graph5.line.bot_nfkb = nanmean(graph5.var_nfkb) - nanstd(graph5.var_nfkb); %creates lower bound
    c = fill([graph5.t,graph5.t(end:-1:1)],[graph5.line.top_nfkb,graph5.line.bot_nfkb(end:-1:1)],colors.sky_blue); %creates filled polygon determined by time axis and upper and lower bound
        set(c, 'EdgeColor', colors.sky_blue);
    alpha(c,0.3)%adds transparency
    hold on
    graph5.line.main_nfkb = plot(graph5.t, nanmean(graph5.var_nfkb)); %adds the mean to the figure
    set(graph5.line.main_nfkb,'Color',colors.sky_blue,'LineWidth',3); %determined properties of the mean line
    hold on
    legend([graph1.line.main_nfkb,graph2.line.main_nfkb,graph3.line.main_nfkb,graph4.line.main_nfkb,graph5.line.main_nfkb],{info1.name,info2.name,info3.name,info4.name,info5.name}, 'Interpreter', 'none')
    
end
if isnumeric(p.Results.id6)    
    graph6.line.top_nfkb = nanmean(graph6.var_nfkb) + nanstd(graph6.var_nfkb); %creates upper bound
    graph6.line.bot_nfkb = nanmean(graph6.var_nfkb) - nanstd(graph6.var_nfkb); %creates lower bound
    c = fill([graph6.t,graph6.t(end:-1:1)],[graph6.line.top_nfkb,graph6.line.bot_nfkb(end:-1:1)],colors.blue); %creates filled polygon determined by time axis and upper and lower bound
        set(c, 'EdgeColor', colors.blue);
    alpha(c,0.3)%adds transparency
    hold on
    graph6.line.main_nfkb = plot(graph6.t, nanmean(graph6.var_nfkb)); %adds the mean to the figure
    set(graph6.line.main_nfkb,'Color',colors.blue,'LineWidth',3); %determined properties of the mean line
    % legend(info1.name,info2.name,info3.name,info4.name,info5.name,info6.name, 'Interpreter', 'none');
    legend([graph1.line.main_nfkb,graph2.line.main_nfkb,graph3.line.main_nfkb,graph4.line.main_nfkb,graph5.line.main_nfkb,graph6.line.main_nfkb],{info1.name,info2.name,info3.name,info4.name,info5.name,info6.name}, 'Interpreter', 'none')
    hold off
    
end
    
    
    % KTR plot
   
    subplot(2,1,2)
    
    graph1.line.top_ktr = nanmean(graph1.var_ktr) + nanstd(graph1.var_ktr); %creates upper bound
    graph1.line.bot_ktr = nanmean(graph1.var_ktr) - nanstd(graph1.var_ktr); %creates lower bound 
    c = fill([graph1.t,graph1.t(end:-1:1)],[graph1.line.top_ktr,graph1.line.bot_ktr(end:-1:1)],colors.dark_red); %creates filled polygon determined by time axis and upper and lower bound
    set(c, 'EdgeColor', colors.dark_red);
    alpha(c,0.4)%adds transparency
    hold on
    graph1.line.main_ktr = plot(graph1.t, nanmean(graph1.var_ktr)); %adds the mean to the figure
    set(graph1.line.main_ktr,'Color', colors.dark_red,'LineWidth',3); %determined properties of the mean line
    hold on
     legend([graph1.line.main_ktr],{info1.name}, 'Interpreter', 'none');
    
    set(gca,'XTick',graph1.opt_ktr.TimeTicks,'YTick',graph1.opt_ktr.MeasurementTicks) %sets axis tick properites
    set(gca,'YTickLabel',graph1.opt_ktr.MeasurementTickLabels,'TickLength',[0.005 0.005]) %sets tick label properties
    ylabel('KTR activation', 'FontSize', 14);
    xlabel('Time (h)','FontSize',14);
    axis([min(graph1.opt_ktr.Times) max(graph1.opt_ktr.Times) graph1.opt_ktr.MeasurementBounds]) %sets axis to determined timescale and determined measurement bounds, derived from maketicks function

if isnumeric(p.Results.id2)    
    graph2.line.top_ktr = nanmean(graph2.var_ktr) + nanstd(graph2.var_ktr); %creates upper bound
    graph2.line.bot_ktr = nanmean(graph2.var_ktr) - nanstd(graph2.var_ktr); %creates lower bound
    c = fill([graph2.t,graph2.t(end:-1:1)],[graph2.line.top_ktr,graph2.line.bot_ktr(end:-1:1)],colors.light_red); %creates filled polygon determined by time axis and upper and lower bound
    set(c, 'EdgeColor', colors.light_red);
    alpha(c,0.5)%adds transparency
    hold on
    graph2.line.main_ktr = plot(graph2.t, nanmean(graph2.var_ktr)); %adds the mean to the figure
    set(graph2.line.main_ktr,'Color',colors.light_red,'LineWidth',3); %determined properties of the mean line
    hold on
    legend([graph1.line.main_ktr,graph2.line.main_ktr],{info1.name,info2.name}, 'Interpreter', 'none')
end

if isnumeric(p.Results.id3)    
    graph3.line.top_ktr = nanmean(graph3.var_ktr) + nanstd(graph3.var_ktr); %creates upper bound
    graph3.line.bot_ktr = nanmean(graph3.var_ktr) - nanstd(graph3.var_ktr); %creates lower bound
    c = fill([graph3.t,graph3.t(end:-1:1)],[graph3.line.top_ktr,graph3.line.bot_ktr(end:-1:1)],colors.orange2); %creates filled polygon determined by time axis and upper and lower bound
    set(c, 'EdgeColor', colors.orange2);
    alpha(c,0.3)%adds transparency
    hold on
    graph3.line.main_ktr = plot(graph3.t, nanmean(graph3.var_ktr)); %adds the mean to the figure
    set(graph3.line.main_ktr,'Color',colors.orange2,'LineWidth',3); %determined properties of the mean line
    hold on
    legend([graph1.line.main_ktr,graph2.line.main_ktr,graph3.line.main_ktr],{info1.name,info2.name,info3.name}, 'Interpreter', 'none')
end
    if isnumeric(p.Results.id4)    
    graph4.line.top_ktr = nanmean(graph4.var_ktr) + nanstd(graph4.var_ktr); %creates upper bound
    graph4.line.bot_ktr = nanmean(graph4.var_ktr) - nanstd(graph4.var_ktr); %creates lower bound
    c = fill([graph4.t,graph4.t(end:-1:1)],[graph4.line.top_ktr,graph4.line.bot_ktr(end:-1:1)],colors.turquoise); %creates filled polygon determined by time axis and upper and lower bound
        set(c, 'EdgeColor', colors.turquoise);
    alpha(c,0.3)%adds transparency
    hold on
    graph4.line.main_ktr = plot(graph4.t, nanmean(graph4.var_ktr)); %adds the mean to the figure
    set(graph4.line.main_ktr,'Color',colors.turquoise,'LineWidth',3); %determined properties of the mean line
    hold on
legend([graph1.line.main_ktr,graph2.line.main_ktr,graph3.line.main_ktr,graph4.line.main_ktr],{info1.name,info2.name,info3.name,info4.name}, 'Interpreter', 'none')
end
if isnumeric(p.Results.id5)    
    graph5.line.top_ktr = nanmean(graph5.var_ktr) + nanstd(graph5.var_ktr); %creates upper bound
    graph5.line.bot_ktr = nanmean(graph5.var_ktr) - nanstd(graph5.var_ktr); %creates lower bound
    c = fill([graph5.t,graph5.t(end:-1:1)],[graph5.line.top_ktr,graph5.line.bot_ktr(end:-1:1)],colors.sky_blue); %creates filled polygon determined by time axis and upper and lower bound
        set(c, 'EdgeColor', colors.sky_blue);
    alpha(c,0.2)%adds transparency
    hold on
    graph5.line.main_ktr = plot(graph5.t, nanmean(graph5.var_ktr)); %adds the mean to the figure
    set(graph5.line.main_ktr,'Color',colors.sky_blue,'LineWidth',3); %determined properties of the mean line
    hold on
    legend([graph1.line.main_ktr,graph2.line.main_ktr,graph3.line.main_ktr,graph4.line.main_ktr,graph5.line.main_ktr],{info1.name,info2.name,info3.name,info4.name,info5.name}, 'Interpreter', 'none')
end
if isnumeric(p.Results.id6)    
    graph6.line.top_ktr = nanmean(graph6.var_ktr) + nanstd(graph6.var_ktr); %creates upper bound
    graph6.line.bot_ktr = nanmean(graph6.var_ktr) - nanstd(graph6.var_ktr); %creates lower bound
    c = fill([graph6.t,graph6.t(end:-1:1)],[graph6.line.top_ktr,graph6.line.bot_ktr(end:-1:1)],colors.blue); %creates filled polygon determined by time axis and upper and lower bound
        set(c, 'EdgeColor', colors.blue);
    alpha(c, 0.2)%adds transparency
    hold on
    graph6.line.main_ktr = plot(graph6.t, nanmean(graph6.var_ktr)); %adds the mean to the figure
    set(graph6.line.main_ktr,'Color',colors.blue,'LineWidth',3); %determined properties of the mean line
    legend([graph1.line.main_ktr,graph2.line.main_ktr,graph3.line.main_ktr,graph4.line.main_ktr,graph5.line.main_ktr,graph6.line.main_ktr],{info1.name,info2.name,info3.name,info4.name,info5.name,info6.name}, 'Interpreter', 'none')
 
    hold off
end
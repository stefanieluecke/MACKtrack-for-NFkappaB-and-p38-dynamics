function [] = see_NFkB_KTR_ratio_HM_sorted_test(id, varargin)

%test visualizing both NFkB and KTR as ratios
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
addParameter(p,'MinLifetime',117, @isnumeric); %allows adjustment of minimum lifetime (?)
addParameter(p,'MinSize',90, valid_conv); %allows adjustment of minimum size (?)
addParameter(p,'TrimFrame',157, @isnumeric);
addParameter(p,'StartThreshNFkB',2, valid_conv); %max allowable starting threshhold to filter out cells with pre-activated NFkB, default is 2
addParameter (p, 'OnThreshNFkB', 0, @isnumeric); %? not used in code?
addParameter (p, 'GraphLimitsNFkB',[-0.05 1],@isnumeric);
addParameter(p,'StartThreshKTR',0.6, valid_conv); %max allowable starting threshhold to filter out cells with pre-activated KTR, default is 0.6
addParameter (p, 'OnThreshKTR', 0, @isnumeric);%? not used in code?
addParameter (p, 'GraphLimitsKTR',[-0.05 0.3],@isnumeric);
addParameter(p, 'SortMetric', 'peakfreq_nfkb');
expectedOrder = {'ascend', 'descend'};
addParameter(p, 'SortOrder', 'descend', @(x)any(validatestring(x, expectedOrder))); 
addParameter(p, 'StartTimePoint', 13, @isnumeric)

parse(p,id, varargin{:})
%%
[metrics,aux, graph, info, measure] = nfkb_ratio_ktr_ratio_metrics_test(id, 'MinLifetime',p.Results.MinLifetime,...
                            'ConvectionShift',p.Results.ConvectionShift, 'OnThreshNFkB',p.Results.OnThreshNFkB,'OnThreshKTR',p.Results.OnThreshKTR,...
                            'MinSize', p.Results.MinSize,'StartThreshNFkB', p.Results.StartThreshNFkB,'StartThreshKTR', p.Results.StartThreshKTR, 'Verbose', ... 
                            p.Results.Verbose, 'GraphLimitsNFkB', p.Results.GraphLimitsNFkB,'GraphLimitsKTR', p.Results.GraphLimitsKTR, 'TrimFrame', p.Results.TrimFrame, 'StartTimePoint', p.Results.StartTimePoint);

SortMetric = p.Results.SortMetric;
[~,graph.order] = sort(metrics.(SortMetric), p.Results.SortOrder);

% Heatmap
figs.a = figure('name','Heatmap_NFkB');
set(figs.a,'Position', [500 7 400 600])
colormapStack(graph.var_nfkb(graph.order,:),graph.celldata(graph.order,:), graph.opt_nfkb); %calls subfunction that makes figure    

%{
%attempt to place both into one graph
pos1 = [300 700 400 600];
subplot('Position', pos1)
colormapStack(graph.var_nfkb(graph.order,:),graph.celldata(graph.order,:), graph.opt_nfkb); %calls subfunction that makes figure    
pos2 = [1300 70 400 600];
subplot('Position', pos2)
colormapStack(graph.var_ktr(graph.order,:),graph.celldata(graph.order,:), graph.opt_ktr); %calls subfunction that makes figure    
%}

figs.b = figure('name','Heatmap_KTR');
set(figs.b,'Position', [500 400 400 600])
colormapStack(graph.var_ktr(graph.order,:),graph.celldata(graph.order,:), graph.opt_ktr); %calls subfunction that makes figure    

figs.e = figure('name', 'Heatmaps_NFkB_KTR');
set(figs.e,'Position', [500 400 1000 600])
colormapStack_double(graph.var_nfkb(graph.order,:), graph.var_ktr(graph.order,:),graph.celldata(graph.order,:), graph.opt_nfkb, graph.opt_ktr);

%{
% try with subplot 
j(1)=subplot(1,2,1);
j(2)=subplot(1,2,2);
% Paste figures on the subplots

copyobj(allchild(get(figs.a,'CurrentAxes')), j(1));
copyobj(allchild(get(figs.b,'CurrentAxes')),j(2));

%copyobj(h1,j(1));
%copyobj(h2,j(2));
%copyobj(allchild(get(h1,'CurrentAxes')),j(1));
%copyobj(allchild(get(h2,'CurrentAxes')),j(2));
%}

%add in line plot of averages and small multiple plots
    colors = setcolors;    

    
% Line plot (mean+/-std) NFkB and KTR next to each other

    graph.line.top = nanmean(graph.var_nfkb) + nanstd(graph.var_nfkb); %creates upper bound
    graph.line.bot = nanmean(graph.var_nfkb) - nanstd(graph.var_nfkb); %creates lower bound

    figs.f = figure('name','MeanTrajectories'); %creates figure and names it
    set(figs.f,'Position', [255   743   1100   500]) %positions figure within window
    subplot(2,1,1);
    fill([graph.t,graph.t(end:-1:1)],[graph.line.top,graph.line.bot(end:-1:1)],colors.light_blue) %creates filled polygon determined by time axis and upper and lower bound
    hold on
    graph.line.main = plot(graph.t, nanmean(graph.var_nfkb)); %adds the mean to the figure
    set(graph.line.main,'Color',[.25 .25 .25],'LineWidth',2); %determined properties of the mean line
    hold off
    set(gca,'XTick',graph.opt_nfkb.TimeTicks,'YTick',graph.opt_nfkb.MeasurementTicks) %sets axis tick properites
    set(gca,'YTickLabel',graph.opt_nfkb.MeasurementTickLabels,'TickLength',[0.005 0.005]) %sets tick label properties
    ylabel('NFkB activation', 'FontSize', 14);
    xlabel('Time (h)','FontSize',14);
    axis([min(graph.opt_nfkb.Times) max(graph.opt_nfkb.Times) graph.opt_nfkb.MeasurementBounds]) %sets axis to determined timescale and determined measurement bounds, derived from maketicks function

    subplot(2,1,2);
    graph.line.top = nanmean(graph.var_ktr) + nanstd(graph.var_ktr); %creates upper bound
    graph.line.bot = nanmean(graph.var_ktr) - nanstd(graph.var_ktr); %creates lower bound
    fill([graph.t,graph.t(end:-1:1)],[graph.line.top,graph.line.bot(end:-1:1)],colors.light_blue) %creates filled polygon determined by time axis and upper and lower bound
    hold on
    graph.line.main = plot(graph.t, nanmean(graph.var_ktr)); %adds the mean to the figure
    set(graph.line.main,'Color',[.25 .25 .25],'LineWidth',2); %determined properties of the mean line
    hold off
    set(gca,'XTick',graph.opt_ktr.TimeTicks,'YTick',graph.opt_ktr.MeasurementTicks) %sets axis tick properites
    set(gca,'YTickLabel',graph.opt_ktr.MeasurementTickLabels,'TickLength',[0.005 0.005]) %sets tick label properties
    ylabel('KTR activation', 'FontSize', 14);
    xlabel('Time (h)','FontSize',14);
    axis([min(graph.opt_ktr.Times) max(graph.opt_ktr.Times) graph.opt_ktr.MeasurementBounds]) %sets axis to determined timescale and determined measurement bounds, derived from maketicks function


% Line plot NFkB(mean+/-std)
graph.line.top = nanmean(graph.var_nfkb) + nanstd(graph.var_nfkb); %creates upper bound
graph.line.bot = nanmean(graph.var_nfkb) - nanstd(graph.var_nfkb); %creates lower bound
figs.c = figure('name','MeanTrajectoryNFkB'); %creates figure and names it
fill([graph.t,graph.t(end:-1:1)],[graph.line.top,graph.line.bot(end:-1:1)],colors.light_blue) %creates filled polygon determined by time axis and upper and lower bound
hold on
graph.line.main = plot(graph.t, nanmean(graph.var_nfkb)); %adds the mean to the figure
set(graph.line.main,'Color',[.25 .25 .25],'LineWidth',2); %determined properties of the mean line
hold off
set(gca,'XTick',graph.opt_nfkb.TimeTicks,'YTick',graph.opt_nfkb.MeasurementTicks) %sets axis tick properites
%set(gca,'YTickLabel',graph.opt_nfkb.MeasurementTickLabels,'TickLength',[0.005 0.005]) %sets tick label properties
set(gca,'YTickLabel',graph.opt_nfkb.MeasurementTickLabels,'TickLength',[0.005 0.005]) %sets tick label properties
ylabel('NFkB activation', 'FontSize', 14);
xlabel('Time (h)','FontSize',14);
set(figs.c,'Position', [555   743   1100   300]) %positions figure within window
axis([min(graph.opt_nfkb.Times) max(graph.opt_nfkb.Times) graph.opt_nfkb.MeasurementBounds]) %sets axis to determined timescale and determined measurement bounds, derived from maketicks function

% Small multiples line graph NFkB
    figs.d  = figure('name','smallmultiplesNFkB');
    set(figs.d,'Position',[500, 350, 876, 1000]);
    graph.smult.order = randperm(size(graph.var_nfkb,1),min([60 size(graph.var_nfkb,1)])); %order the plots randomly (vector of max 60 or number of cells random integers between 1 and the number of cells)
    graph.smult.h = tight_subplot(15,4); %creates 15x4 small subplots using function
    xpos = max(graph.opt_nfkb.Times)-0.48*(max(graph.opt_nfkb.Times)-min(graph.opt_nfkb.Times)); %determines the positioning of text within each small graph
    ypos =  max(graph.opt_nfkb.MeasurementBounds) - 0.26*diff(graph.opt_nfkb.MeasurementBounds); %determines the positioning of text within each small graph
    for i =1:length(graph.smult.order)
        %plot(graph.smult.h(i),graph.t,graph.var(graph.smult.order(i),:),'Color',colors.grays{3}, 'LineWidth',2)
        plot(graph.smult.h(i),graph.t,graph.var_nfkb(graph.smult.order(i),:),'Color',colors.red, 'LineWidth',2)
        set(graph.smult.h(i),'XLim',[min(graph.opt_nfkb.Times)-(range(graph.opt_nfkb.Times)*0.02) max(graph.opt_nfkb.Times)],...
            'YLim',graph.opt_nfkb.MeasurementBounds,'XTickLabel',{},'YTickLabel',{})
        text(xpos,ypos,['XY ',num2str(graph.celldata(graph.smult.order(i),1)),...
            ', cell ',num2str(graph.celldata(graph.smult.order(i),2))],'Parent',graph.smult.h(i))
            
    end
% Line plot KTR (mean+/-std)
    graph.line.top = nanmean(graph.var_ktr) + nanstd(graph.var_ktr); %creates upper bound
    graph.line.bot = nanmean(graph.var_ktr) - nanstd(graph.var_ktr); %creates lower bound
    figs.c = figure('name','MeanTrajectoryKTR'); %creates figure and names it
    fill([graph.t,graph.t(end:-1:1)],[graph.line.top,graph.line.bot(end:-1:1)],colors.light_blue) %creates filled polygon determined by time axis and upper and lower bound
    hold on
    graph.line.main = plot(graph.t, nanmean(graph.var_ktr)); %adds the mean to the figure
    set(graph.line.main,'Color',[.25 .25 .25],'LineWidth',2); %determined properties of the mean line
    hold off
    set(gca,'XTick',graph.opt_ktr.TimeTicks,'YTick',graph.opt_ktr.MeasurementTicks) %sets axis tick properites
    set(gca,'YTickLabel',graph.opt_ktr.MeasurementTickLabels,'TickLength',[0.005 0.005]) %sets tick label properties
    set(figs.c,'Position', [555   743   1100   300]) %positions figure within window
    ylabel('KTR activation', 'FontSize', 14);
    xlabel('Time (h)','FontSize',14);
    axis([min(graph.opt_ktr.Times) max(graph.opt_ktr.Times) graph.opt_ktr.MeasurementBounds]) %sets axis to determined timescale and determined measurement bounds, derived from maketicks function

% Small multiples line graph KTR
    figs.d  = figure('name','smallmultiplesKTR');
    set(figs.d,'Position',[500, 350, 876, 1000]);
    graph.smult.order = randperm(size(graph.var_ktr,1),min([60 size(graph.var_ktr,1)])); %order the plots randomly (vector of max 60 or number of cells random integers between 1 and the number of cells)
    graph.smult.h = tight_subplot(15,4); %creates 15x4 small subplots using function
    xpos = max(graph.opt_ktr.Times)-0.48*(max(graph.opt_ktr.Times)-min(graph.opt_ktr.Times)); %determines the positioning of text within each small graph
    ypos =  max(graph.opt_ktr.MeasurementBounds) - 0.26*diff(graph.opt_ktr.MeasurementBounds); %determines the positioning of text within each small graph
    for i =1:length(graph.smult.order)
        %plot(graph.smult.h(i),graph.t,graph.var(graph.smult.order(i),:),'Color',colors.grays{3}, 'LineWidth',2)
        plot(graph.smult.h(i),graph.t,graph.var_ktr(graph.smult.order(i),:),'Color',colors.red, 'LineWidth',2)
        set(graph.smult.h(i),'XLim',[min(graph.opt_ktr.Times)-(range(graph.opt_ktr.Times)*0.02) max(graph.opt_ktr.Times)],...
            'YLim',graph.opt_ktr.MeasurementBounds,'XTickLabel',{},'YTickLabel',{})
        text(xpos,ypos,['XY ',num2str(graph.celldata(graph.smult.order(i),1)),...
            ', cell ',num2str(graph.celldata(graph.smult.order(i),2))],'Parent',graph.smult.h(i))
            
    end
end

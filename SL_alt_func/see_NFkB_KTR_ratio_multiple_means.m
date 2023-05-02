function [ID] = see_NFkB_KTR_ratio_multiple_means(IDs, varargin)

% 
%Enter highest responder first
%% INPUT PARSING
% Create input parser object, add required params from function input
p = inputParser;
% Required: ID input
valid_id = @(x) assert((isnumeric(x)&&length(x)==1)||exist(x,'file'),...
    'ID input must be spreadsheet ID or full file path');
addRequired(p,'IDs');

expectedFlags = {'on','off'};
addParameter(p,'Verbose','off', @(x) any(validatestring(x,expectedFlags)));%checks whether optional name-value argument matches on or off %checks if x matches expectedFlags
valid_conv = @(x) assert(isnumeric(x)&&(x>=0)&&(length(x)==1),...
    'Parameter must be single integer >= 0'); %checks whether parameters below are single integers
addParameter(p,'ConvectionShift',1, valid_conv); %allows adjustment of convection shift (?)
addParameter(p,'MinLifetime',109, @isnumeric); %allows adjustment of minimum lifetime (?)
addParameter(p,'MinSize',90, valid_conv); %allows adjustment of minimum size (?)
addParameter(p,'TrimFrame',157, @isnumeric);
addParameter(p,'StartThreshNFkB',14, valid_conv); %max allowable starting threshhold to filter out cells with pre-activated NFkB, default is 2
addParameter (p, 'OnThreshNFkB', 3, @isnumeric); %? not used in code?
addParameter (p, 'GraphLimitsNFkB',[-0.9 7],@isnumeric);
addParameter(p,'StartThreshKTR',0.9, valid_conv); %max allowable starting threshhold to filter out cells with pre-activated KTR, default is 0.6
addParameter (p, 'OnThreshKTR', 3, @isnumeric);%? not used in code?
addParameter (p, 'GraphLimitsKTR',[-0.02 0.35],@isnumeric);
addParameter(p, 'StimulationTimePoint', 13, @isnumeric)
addParameter(p, 'FramesPerHour', 12, @isnumeric)
addParameter(p,'NFkBBaselineDeduction', 'on', @(x) any(validatestring(x,expectedFlags))) %option to turn off NFkB baseline deduction
addParameter(p, 'NFkBBackgroundAdjustment', 'on',@(x) any(validatestring(x,expectedFlags))) %option to turn off NFkB fluorescence distribution adjustment
addParameter(p,'NFkBBaselineAdjustment', 'on', @(x) any(validatestring(x,expectedFlags))) %option to turn off adjusment of NFkB trajectories with correction factor for fluorescence drop derived from Mock experiments
addParameter(p,'KTRBaselineDeduction', 'on', @(x) any(validatestring(x,expectedFlags))) %option to turn off NFkB baseline deduction
addParameter(p,'KTRBaselineAdjustment', 'on', @(x) any(validatestring(x,expectedFlags))) %option to turn off adjusment of KTR trajectories with correction factor for fluorescence drop derived from Mock experiments



expectedFilters = {'none','nfkb', 'ktr', 'both', 'respective'};
addParameter(p, 'FilterResponders','none', @(x) any(validatestring(x,expectedFilters)));%filter out non-responders or not
addParameter(p, 'IncludeKTR', 'on',@(x)any(validatestring(x, expectedFlags)));

parse(p,IDs, varargin{:})

colors = setcolors;    

n = numel(IDs);
ID(n).metrics = [];
ID(n).graph = [];
ID(n).info = [];
ID(n).aux = [];
ID(n).measure = [];

%%
for i= 1:n
[ID(i).metrics,ID(i).aux, ID(i).graph, ID(i).info, ID(i).measure] = nfkb_ktr_ratio_metrics(IDs(i), 'MinLifetime',p.Results.MinLifetime,...
                             'OnThreshNFkB',p.Results.OnThreshNFkB,'OnThreshKTR',p.Results.OnThreshKTR,...
                            'MinSize', p.Results.MinSize, 'Verbose', ... 
                            p.Results.Verbose, 'GraphLimitsNFkB', p.Results.GraphLimitsNFkB,'GraphLimitsKTR', p.Results.GraphLimitsKTR, 'TrimFrame', p.Results.TrimFrame,...
                            'StimulationTimePoint', p.Results.StimulationTimePoint,'FramesPerHour', p.Results.FramesPerHour, 'NFkBBaselineDeduction', p.Results.NFkBBaselineDeduction, 'NFkBBackgroundAdjustment',p.Results.NFkBBackgroundAdjustment,...
                            'NFkBBaselineAdjustment', p.Results.NFkBBaselineAdjustment,'KTRBaselineDeduction', p.Results.KTRBaselineDeduction,'KTRBaselineAdjustment', p.Results.KTRBaselineAdjustment,'IncludeKTR',p.Results.IncludeKTR);
end


switch p.Results.FilterResponders 
    case 'nfkb'
        for i = 1:n
           ID(i).graph.var_nfkb     = ID(i).graph.var_nfkb(ID(i).metrics.responder_index_nfkb == 1,:);
           ID(i).graph.var_ktr      = ID(i).graph.var_ktr(ID(i).metrics.responder_index_nfkb == 1,:);
                   end
    case 'ktr'
        for i = 1:n
           ID(i).graph.var_nfkb     = ID(i).graph.var_nfkb(ID(i).metrics.responder_index_ktr == 1,:);
           ID(i).graph.var_ktr      = ID(i).graph.var_ktr(ID(i).metrics.responder_index_ktr == 1,:);
        end
    case 'both'
        for i = 1:n
           ID(i).graph.var_nfkb     = ID(i).graph.var_nfkb((ID(i).metrics.responder_index_nfkb == 1     & ID(i).metrics.responder_index_ktr == 1),:);
           ID(i).graph.var_ktr      = ID(i).graph.var_ktr((ID(i).metrics.responder_index_nfkb == 1      & ID(i).metrics.responder_index_ktr == 1),:);
        end
    case 'respective'
        for i = 1:n
           ID(i).graph.var_nfkb     = ID(i).graph.var_nfkb(ID(i).metrics.responder_index_nfkb == 1,:);
           ID(i).graph.var_ktr      = ID(i).graph.var_ktr(ID(i).metrics.responder_index_ktr == 1,:);
        end
end

%%
% Line plot NFkB(mean+/-std)
if strcmpi(p.Results.IncludeKTR,'on')

    %Line plot NFkB and KTR (mean+/-std)
    figs.a = figure('name','OverlayMeanTrajectoriesNFkBandKTR'); %creates figure and names it
    set(figs.a,'Position', [555   743   750   700]) %positions figure within window
    %NFkB plot
    subplot(2,1,1)   
        text_for_legend = [];
        lines_for_legend = [];
        for i = n:-1:1 %this insures the lower ID numbers (usually lower dose are plotted on top of the higher doses and not hidden
            ID(i).graph.line.top_nfkb = nanmean(ID(i).graph.var_nfkb) + nanstd(ID(i).graph.var_nfkb); %creates upper bound
            ID(i).graph.line.bot_nfkb = nanmean(ID(i).graph.var_nfkb) - nanstd(ID(i).graph.var_nfkb); %creates lower bound
            c(i) = fill([ID(i).graph.t,ID(i).graph.t(end:-1:1)],[ID(i).graph.line.top_nfkb,ID(i).graph.line.bot_nfkb(end:-1:1)], colors.doses{mod(i-1,length(colors.doses))+1}); %creates filled polygon determined by time axis and upper and lower bound
            set(c(i), 'EdgeColor',colors.doses{mod(i-1,length(colors.doses))+1});
            alpha(c(i),0.1)%adds transparency
            hold on
            ID(i).graph.line.main_nfkb = plot( ID(i).graph.t, nanmean( ID(i).graph.var_nfkb)); %adds the mean to the figure
            set( ID(i).graph.line.main_nfkb,'Color', colors.doses{mod(i-1,length(colors.doses))+1},'LineWidth',3); %determined properties of the mean line
            text_for_legend = [{ID(i).info.name},text_for_legend];
            lines_for_legend = [ ID(i).graph.line.main_nfkb,lines_for_legend];
        end

        set(gca,'XTick',ID(n).graph.opt_nfkb.TimeTicks,'YTick',ID(n).graph.opt_nfkb.MeasurementTicks) %sets axis tick properites
        set(gca,'YTickLabel',ID(n).graph.opt_nfkb.MeasurementTickLabels,'TickLength',[0.005 0.005]) %sets tick label properties
        ylabel('NFkB activation', 'FontSize', 14);
        xlabel('Time (h)','FontSize',14);
        axis([min(ID(n).graph.opt_nfkb.Times) max(ID(n).graph.opt_nfkb.Times) ID(n).graph.opt_nfkb.MeasurementBounds]) %sets axis to determined timescale and determined measurement bounds, derived from maketicks function
        legend(lines_for_legend,cellstr(text_for_legend), 'Interpreter', 'none')  
        hold off

    %KTR plot
    subplot(2,1,2)
        text_for_legend = [];
        lines_for_legend = [];
            for i = n:-1:1 %this insures the lower ID numbers (usually lower dose are plotted on top of the higher doses and not hidden
                ID(i).graph.line.top_ktr = nanmean(ID(i).graph.var_ktr) + nanstd(ID(i).graph.var_ktr); %creates upper bound
                ID(i).graph.line.bot_ktr = nanmean(ID(i).graph.var_ktr) - nanstd(ID(i).graph.var_ktr); %creates lower bound
                b(i) = fill([ID(i).graph.t,ID(i).graph.t(end:-1:1)],[ID(i).graph.line.top_ktr,ID(i).graph.line.bot_ktr(end:-1:1)], colors.doses{mod(i-1,length(colors.doses))+1}); %creates filled polygon determined by time axis and upper and lower bound
                set(b(i), 'EdgeColor',colors.doses{mod(i-1,length(colors.doses))+1});
                alpha(b(i),0.1)%adds transparency
                hold on
                ID(i).graph.line.main_ktr = plot( ID(i).graph.t, nanmean( ID(i).graph.var_ktr)); %adds the mean to the figure
                set( ID(i).graph.line.main_ktr,'Color', colors.doses{mod(i-1,length(colors.doses))+1},'LineWidth',3); %determined properties of the mean line
                text_for_legend = [{ID(i).info.name},text_for_legend];
                lines_for_legend = [ ID(i).graph.line.main_ktr,lines_for_legend];
            end

        set(gca,'XTick',ID(n).graph.opt_ktr.TimeTicks,'YTick',ID(n).graph.opt_ktr.MeasurementTicks) %sets axis tick properites
        set(gca,'YTickLabel',ID(n).graph.opt_ktr.MeasurementTickLabels,'TickLength',[0.005 0.005]) %sets tick label properties
        ylabel('KTR activation', 'FontSize', 14);
        xlabel('Time (h)','FontSize',14);
        axis([min(ID(n).graph.opt_ktr.Times) max(ID(n).graph.opt_ktr.Times) ID(n).graph.opt_ktr.MeasurementBounds]) %sets axis to determined timescale and determined measurement bounds, derived from maketicks function
        legend(lines_for_legend,cellstr(text_for_legend), 'Interpreter', 'none')  
        hold off

    %
else
   % Only Line plot NFkB(mean+/-std)

    figs.c = figure('name','OverlayMeanTrajectoriesNFkB'); %creates figure and names it
    set(figs.c,'Position', [555   743   750   350]) %positions figure within window
    text_for_legend = [];
    lines_for_legend = [];
        for i = n:-1:1 %this insures the lower ID numbers (usually lower dose are plotted on top of the higher doses and not hidden
            ID(i).graph.line.top_nfkb = nanmean(ID(i).graph.var_nfkb) + nanstd(ID(i).graph.var_nfkb); %creates upper bound
            ID(i).graph.line.bot_nfkb = nanmean(ID(i).graph.var_nfkb) - nanstd(ID(i).graph.var_nfkb); %creates lower bound
            c(i) = fill([ID(i).graph.t,ID(i).graph.t(end:-1:1)],[ID(i).graph.line.top_nfkb,ID(i).graph.line.bot_nfkb(end:-1:1)], colors.doses{mod(i-1,length(colors.doses))+1}); %creates filled polygon determined by time axis and upper and lower bound
            set(c(i), 'EdgeColor',colors.doses{mod(i-1,length(colors.doses))+1});
            alpha(c(i),0.1)%adds transparency
            hold on
            ID(i).graph.line.main_nfkb = plot( ID(i).graph.t, nanmean( ID(i).graph.var_nfkb)); %adds the mean to the figure
            set( ID(i).graph.line.main_nfkb,'Color', colors.doses{mod(i-1,length(colors.doses))+1},'LineWidth',3); %determined properties of the mean line
            text_for_legend = [{ID(i).info.name},text_for_legend];
            lines_for_legend = [ ID(i).graph.line.main_nfkb,lines_for_legend];
        end

    set(gca,'XTick',ID(n).graph.opt_nfkb.TimeTicks,'YTick',ID(n).graph.opt_nfkb.MeasurementTicks) %sets axis tick properites
    set(gca,'YTickLabel',ID(n).graph.opt_nfkb.MeasurementTickLabels,'TickLength',[0.005 0.005]) %sets tick label properties
    ylabel('NFkB activation', 'FontSize', 14);
    xlabel('Time (h)','FontSize',14);
    axis([min(ID(n).graph.opt_nfkb.Times) max(ID(n).graph.opt_nfkb.Times) ID(n).graph.opt_nfkb.MeasurementBounds]) %sets axis to determined timescale and determined measurement bounds, derived from maketicks function
    legend(lines_for_legend,cellstr(text_for_legend), 'Interpreter', 'none')  
    hold off
end

%{
%% Line plot KTR (mean+/-std)
    figs.b = figure('name','OverlayMeanTrajectoriesKTR'); %creates figure and names it
    set(figs.b,'Position', [555   743   1145   500]) %positions figure within window
    text_for_legend = [];
    lines_for_legend = [];
        for i = n:-1:1 %this insures the lower ID numbers (usually lower dose are plotted on top of the higher doses and not hidden
            ID(i).graph.line.top_ktr = nanmean(ID(i).graph.var_ktr) + nanstd(ID(i).graph.var_ktr); %creates upper bound
            ID(i).graph.line.bot_ktr = nanmean(ID(i).graph.var_ktr) - nanstd(ID(i).graph.var_ktr); %creates lower bound
            b(i) = fill([ID(i).graph.t,ID(i).graph.t(end:-1:1)],[ID(i).graph.line.top_ktr,ID(i).graph.line.bot_ktr(end:-1:1)], colors.doses{mod(i-1,length(colors.doses))+1}); %creates filled polygon determined by time axis and upper and lower bound
            set(b(i), 'EdgeColor',colors.doses{mod(i-1,length(colors.doses))+1});
            alpha(b(i),0.2)%adds transparency
            hold on
            ID(i).graph.line.main_ktr = plot( ID(i).graph.t, nanmean( ID(i).graph.var_ktr)); %adds the mean to the figure
            set( ID(i).graph.line.main_ktr,'Color', colors.doses{mod(i-1,length(colors.doses))+1},'LineWidth',3); %determined properties of the mean line
            text_for_legend = [{ID(i).info.name},text_for_legend];
            lines_for_legend = [ ID(i).graph.line.main_ktr,lines_for_legend];
        end

    set(gca,'XTick',ID(n).graph.opt_ktr.TimeTicks,'YTick',ID(n).graph.opt_ktr.MeasurementTicks) %sets axis tick properites
    set(gca,'YTickLabel',ID(n).graph.opt_ktr.MeasurementTickLabels,'TickLength',[0.005 0.005]) %sets tick label properties
    ylabel('KTR activation', 'FontSize', 14);
    xlabel('Time (h)','FontSize',14);
    axis([min(ID(n).graph.opt_ktr.Times) max(ID(n).graph.opt_ktr.Times) ID(n).graph.opt_ktr.MeasurementBounds]) %sets axis to determined timescale and determined measurement bounds, derived from maketicks function
    legend(lines_for_legend,cellstr(text_for_legend), 'Interpreter', 'none')  
    hold off

%


%}
function [ID, violin] = graphFeatures(IDs, varargin)

%20200625 graphs features for Violin Plots, %responders, %oscillators for NFkB and KTR
%also provides features with option for filtering based on responder status


%for any number of experiment IDs
p = inputParser;
addRequired(p,'IDs'); %vector containing IDs to be plotted
%parameters to be passed to metric function
expectedFlags = {'on','off'};
valid_conv = @(x) assert(isnumeric(x)&&(x>=0)&&(length(x)==1),...
    'Parameter must be single integer >= 0'); %checks whether parameters below are single integers
addParameter(p,'Verbose','off', @(x) any(validatestring(x,expectedFlags)));%checks whether optional name-value argument matches on or off %checks if x matches expectedFlags
addParameter(p,'MinLifetime',109, @isnumeric); %allows adjustment of minimum lifetime (?)
addParameter(p,'MinSize',90, valid_conv); %allows adjustment of minimum size (?)
addParameter(p,'TrimFrame',157, @isnumeric);
addParameter(p,'StartThreshNFkB',14, valid_conv);%max allowable starting threshhold (before baseline deduction)to filter out cells with pre-activated NFkB
addParameter (p, 'OnThreshNFkB', 3, @isnumeric); %sigma threshold for determining responders
addParameter (p, 'GraphLimitsNFkB',[-0.25 7],@isnumeric);
addParameter(p,'StartThreshKTR',0.9, valid_conv); %max allowable starting threshhold to filter out cells with pre-activated KTR, default is 0.6
addParameter (p, 'OnThreshKTR', 3, @isnumeric); %sigma threshold for determining responders
addParameter (p, 'GraphLimitsKTR',[-0.02,0.35],@isnumeric);
addParameter(p, 'StimulationTimePoint', 13, @isnumeric);
addParameter(p, 'FramesPerHour', 12, @isnumeric);
addParameter(p,'NFkBBaselineAdjustment', 'on', @(x) any(validatestring(x,expectedFlags))) %option to turn off adjusment of NFkB trajectories with correction factor for fluorescence drop derived from Mock experiments
addParameter(p, 'IncludeKTR', 'on',@(x)any(validatestring(x, expectedFlags)));

addParameter(p, 'MetricFilter', 'off',@(x)any(validatestring(x, expectedFlags)));
addParameter(p, 'FilterMetric', 'max_amplitude_4h_nfkb1');
addParameter (p, 'MetricThresh', 1, @isnumeric);

addParameter (p, 'RestrictCellNumber', 'off',@(x)any(validatestring(x, expectedFlags)));
addParameter (p, 'CellNumber', []);

%parameter to access metrics to be graphed in violin plots
addParameter(p, 'FeatureListFile', 'C:\Users\stlue\OneDrive\PostDoc UCLA\1 Post Doc UCLA\Matlab analysis\MACKtrack_SL\Metrics_SL\FeatureList.xlsx') %provide file path for Excel table with list of feature to be computed
%addParameter(p, 'FeatureListFile', 'D:\OneDrive\PostDoc UCLA\1 Post Doc UCLA\Matlab analysis\MACKtrack_SL\Metrics_SL\FeatureList.xlsx') %provide file path for Excel table with list of feature to be computed
addParameter(p, 'XLabel', {})
addParameter(p, 'Color', {[1,0,0], [0,1,0]})

expectedFilters = {'none','nfkb', 'ktr', 'both', 'respective'};
addParameter(p, 'FilterResponders','none', @(x) any(validatestring(x,expectedFilters)));%filter out non-responders or not

parse(p,IDs, varargin{:})

FilterMetric = p.Results.FilterMetric;

colors = setcolors; 

ViolMetTableNFkB = readtable(p.Results.FeatureListFile, 'Sheet', 'ViolFeatTableNFkB');
viol_met_nfkb = table2cell(ViolMetTableNFkB(:,1)); %list of metrics/features to be plotted
viol_met_index_nfkb = table2cell(ViolMetTableNFkB(:,2)); %index for metrics with multiple columns (eg duration, etc.)
viol_met_units_nfkb = table2cell(ViolMetTableNFkB(:,3));
ViolMetTableKTR= readtable(p.Results.FeatureListFile, 'Sheet', 'ViolFeatTableKTR');
viol_met_ktr = table2cell(ViolMetTableKTR(:,1)); %list of metrics/features to be plotted
viol_met_index_ktr = table2cell(ViolMetTableKTR(:,2)); %index for metrics with multiple columns (eg duration, etc.)
viol_met_units_ktr = table2cell(ViolMetTableKTR(:,3));

%viol_met_nfkb = {'off_times_nfkb','max_amplitude_nfkb','max_integral_nfkb','pk1_amp_nfkb', 'pk2_amp_nfkb','peakfreq_nfkb','max_derivative_nfkb', 'min_derivative_nfkb', 'pk1_time_nfkb', 'pk2_time_nfkb'};
%viol_met_ktr = {'off_times_ktr','max_amplitude_ktr','max_integral_ktr','pk1_amp_ktr', 'pk2_amp_ktr','peakfreq_ktr','max_derivative_ktr', 'min_derivative_ktr', 'pk1_time_ktr', 'pk2_time_ktr'};

n = numel(IDs);
ID(n).metrics = [];
ID(n).graph = [];
ID(n).info = [];
ID(n).features = [];


% Metrics for NFkB
%run metrics function of desired exp ID
for i= 1:n
    [ID(i).metrics,~,ID(i).graph,ID(i).info,~] = nfkb_ktr_ratio_metrics(IDs(i), 'MinLifetime',p.Results.MinLifetime,...
                            'OnThreshNFkB',p.Results.OnThreshNFkB,'OnThreshKTR',p.Results.OnThreshKTR,...
                            'MinSize', p.Results.MinSize,'StartThreshNFkB', p.Results.StartThreshNFkB,'StartThreshKTR', p.Results.StartThreshKTR, 'Verbose', ... 
                            p.Results.Verbose, 'TrimFrame', p.Results.TrimFrame, 'StimulationTimePoint', p.Results.StimulationTimePoint, 'FramesPerHour', p.Results.FramesPerHour,'IncludeKTR',p.Results.IncludeKTR,'NFkBBaselineAdjustment', p.Results.NFkBBaselineAdjustment);
    if strcmpi(p.Results.IncludeKTR,'on')

        ID(i).features = computeFeatures(IDs(i), 'metrics', ID(i).metrics, 'FeatureListTable', [ViolMetTableNFkB;ViolMetTableKTR],... 
                                    'MinLifetime',p.Results.MinLifetime, 'OnThreshNFkB',p.Results.OnThreshNFkB,'OnThreshKTR',p.Results.OnThreshKTR,...
                                    'MinSize', p.Results.MinSize,'StartThreshNFkB', p.Results.StartThreshNFkB,'StartThreshKTR', p.Results.StartThreshKTR, 'Verbose', ... 
                                    p.Results.Verbose, 'TrimFrame', p.Results.TrimFrame, ...
                                    'StimulationTimePoint', p.Results.StimulationTimePoint, 'FramesPerHour', p.Results.FramesPerHour,'IncludeKTR',p.Results.IncludeKTR);
        else
            ID(i).features = computeFeatures(IDs(i), 'metrics', ID(i).metrics, 'FeatureListTable', [ViolMetTableNFkB],... 
                                    'MinLifetime',p.Results.MinLifetime, 'OnThreshNFkB',p.Results.OnThreshNFkB,'OnThreshKTR',p.Results.OnThreshKTR,...
                                    'MinSize', p.Results.MinSize,'StartThreshNFkB', p.Results.StartThreshNFkB,'StartThreshKTR', p.Results.StartThreshKTR, 'Verbose', ... 
                                    p.Results.Verbose, 'TrimFrame', p.Results.TrimFrame, ...
                                    'StimulationTimePoint', p.Results.StimulationTimePoint, 'FramesPerHour', p.Results.FramesPerHour,'IncludeKTR',p.Results.IncludeKTR);
    end

end

%% Responder Filtering
switch p.Results.FilterResponders 
    case 'none'
       for i = 1:n
        for k = 1:numel(viol_met_nfkb)
            if i==1 
            violin.([viol_met_nfkb{k},num2str(viol_met_index_nfkb{k})]) = {ID(i).features.(viol_met_nfkb{k})(:, viol_met_index_nfkb{k})};
            else
            violin.([viol_met_nfkb{k},num2str(viol_met_index_nfkb{k})]) = {violin.([viol_met_nfkb{k},num2str(viol_met_index_nfkb{k})]){:}, ID(i).features.(viol_met_nfkb{k})(:, viol_met_index_nfkb{k})};
            end
        end
       end
    case {'nfkb', 'respective'}
       for i = 1:n
        for k = 1:numel(viol_met_nfkb)
            if i==1 
            violin.([viol_met_nfkb{k},num2str(viol_met_index_nfkb{k})]) = {ID(i).features.(viol_met_nfkb{k})(ID(i).metrics.responder_index_nfkb == 1, viol_met_index_nfkb{k})};
            else
            violin.([viol_met_nfkb{k},num2str(viol_met_index_nfkb{k})]) = {violin.([viol_met_nfkb{k},num2str(viol_met_index_nfkb{k})]){:}, ID(i).features.(viol_met_nfkb{k})(ID(i).metrics.responder_index_nfkb == 1,viol_met_index_nfkb{k})};
            end
        end
       end
    case 'ktr'
       for i = 1:n
        for k = 1:numel(viol_met_nfkb)
            if i==1 
            violin.([viol_met_nfkb{k},num2str(viol_met_index_nfkb{k})]) = {ID(i).features.(viol_met_nfkb{k})(ID(i).metrics.responder_index_ktr == 1, viol_met_index_nfkb{k})};
            else
            violin.([viol_met_nfkb{k},num2str(viol_met_index_nfkb{k})]) = {violin.([viol_met_nfkb{k},num2str(viol_met_index_nfkb{k})]){:}, ID(i).features.(viol_met_nfkb{k})(ID(i).metrics.responder_index_ktr == 1, viol_met_index_nfkb{k})};
            end
        end
       end
    case 'both'
       for i = 1:n
        for k = 1:numel(viol_met_nfkb)
            if i==1 
            violin.([viol_met_nfkb{k},num2str(viol_met_index_nfkb{k})]) = {ID(i).features.(viol_met_nfkb{k})(ID(i).metrics.responder_index_nfkb == 1 & ID(i).metrics.responder_index_ktr == 1, viol_met_index_nfkb{k})};
            else
            violin.([viol_met_nfkb{k},num2str(viol_met_index_nfkb{k})]) = {violin.([viol_met_nfkb{k},num2str(viol_met_index_nfkb{k})]){:}, ID(i).features.(viol_met_nfkb{k})(ID(i).metrics.responder_index_nfkb == 1 & ID(i).metrics.responder_index_ktr == 1, viol_met_index_nfkb{k})};
            end
        end
       end   
end

%% Filter for a metric threshold if desired (rarely)
if strcmpi(p.Results.MetricFilter,'on')
        for i = 1:n
            for k = 1:numel(viol_met_nfkb)
                    violin.([viol_met_nfkb{k},num2str(viol_met_index_nfkb{k})]){:,i}     = violin.([viol_met_nfkb{k},num2str(viol_met_index_nfkb{k})]){i}(violin.(FilterMetric){:,i} >= p.Results.MetricThresh); 
            end
        end
end

%% Show same number of cells (if desired)

if strcmpi(p.Results.RestrictCellNumber,'on')
     CellNumber = p.Results.CellNumber;
    if isempty(p.Results.CellNumber) 
        for i = 1:n
            CellNumber = [CellNumber, size(violin.([viol_met_nfkb{k},num2str(viol_met_index_nfkb{k})]){i},1)];
        end
        CellNumber = min(CellNumber);
    end

    for i = 1:n
           index = randperm(size(violin.([viol_met_nfkb{k},num2str(viol_met_index_nfkb{k})]){i},1),CellNumber);
           for k = 1:numel(viol_met_nfkb)
                    violin.([viol_met_nfkb{k},num2str(viol_met_index_nfkb{k})]){:,i}     = violin.([viol_met_nfkb{k},num2str(viol_met_index_nfkb{k})]){i}(index); 
           end
     end
end

%% Calculate statistical difference between features % NOTE only works for comparison of 2 cond
for j = 1:numel(viol_met_nfkb)
    if ~iscategorical (ID(i).features.(viol_met_nfkb{j}))
        [~, ks_pval.([viol_met_nfkb{j},num2str(viol_met_index_nfkb{j})]), ~] = kstest2(violin.([viol_met_nfkb{j},num2str(viol_met_index_nfkb{j})]){1}, violin.([viol_met_nfkb{j},num2str(viol_met_index_nfkb{j})]){2}); %continuous 2-sample ks test
    end
end

%% Make violin plots

figure
if strcmpi(p.Results.IncludeKTR,'on')

    if numel(viol_met_nfkb)<7 
        tiledlayout(2,round(numel(viol_met_nfkb))) %looks best when using even number of metrics/features
    else
        tiledlayout(4,round(numel(viol_met_nfkb)/2)) 
    end 
else
    if numel(viol_met_nfkb)<6 
    tiledlayout(1,round(numel(viol_met_nfkb))) %looks best when using even number of metrics/features
    elseif numel(viol_met_nfkb)<12 
    tiledlayout(2,round(numel(viol_met_nfkb)/2)) %looks best when using even number of metrics/features  
    else
    tiledlayout(3,round(numel(viol_met_nfkb)/3)) 
    end 
end

violin_spacing = 1:n;
for k = 1:numel(viol_met_nfkb)
     if ~iscategorical (ID(i).features.(viol_met_nfkb{k}))
     %   violin_mack(violin.(viol_met_nfkb{k}), violin_spacing, 'Area', 0.05, 'YLim', [-5 15])
        axes.ax(k) = nexttile;
        violin_mack_update(violin.([viol_met_nfkb{k},num2str(viol_met_index_nfkb{k})]),violin_spacing,'Axes', axes.ax(k), 'Area', 0.15,'XSpace', 0.4, 'BinScale', 1,'Smoothing', 'on', 'Connect', 'off', 'MarkerSize', 6, 'ShowBins', 'off','Mode', 'off', 'Color', p.Results.Color);
    %   violin_kernel(violin.(viol_met_nfkb{k}),violin_spacing,'Axes', axes.ax(k), 'Area', 0.04,'XSpace', 0.1, 'BinScale', 1,'Smoothing', 'on', 'Connect', 'on', 'MarkerSize', 7, 'ShowBins', 'off');
        set(axes.ax(k),'XTick',1:n, 'XTickLabel',p.Results.XLabel,'Box','on','GridLineStyle', ':','LineWidth',1,...  
       'YTickLabelMode', 'auto', 'Color', ones(1,3)*1,'XGrid', 'on','GridColor', 'w', 'XTickLabelRotation', 45);
        title([viol_met_nfkb{k},' ',num2str(viol_met_index_nfkb{k})], 'Interpreter', 'none')
        ylabel(viol_met_units_nfkb{k})
        if max(ylim)>=0
            text(1.05, max(ylim)*0.95, ['p=' num2str(round(ks_pval.([viol_met_nfkb{k},num2str(viol_met_index_nfkb{k})]),2, 'significant'))])
        elseif max(ylim)< 0
            text(1.05, min(ylim)*0.9, ['p=' num2str(round(ks_pval.([viol_met_nfkb{k},num2str(viol_met_index_nfkb{k})]),2, 'significant'))])
        end
     end
end

%% KTR metrics

if strcmpi(p.Results.IncludeKTR,'on')


    %figure
    %tiledlayout(2,round(numel(viol_met)/2))

    switch p.Results.FilterResponders 
        case 'none'
           for i = 1:n
            for k = 1:numel(viol_met_ktr)
                if i==1 
                violin.(viol_met_ktr{k}) = {ID(i).features.(viol_met_ktr{k})(:, viol_met_index_ktr{k})};
                else
                violin.(viol_met_ktr{k}) = [violin.(viol_met_ktr{k}), ID(i).features.(viol_met_ktr{k})(:, viol_met_index_ktr{k})];
                end
            end
           end
        case 'nfkb'
           for i = 1:n
            for k = 1:numel(viol_met_ktr)
                if i==1 
                violin.(viol_met_ktr{k}) = {ID(i).features.(viol_met_ktr{k})(ID(i).metrics.responder_index_nfkb == 1,viol_met_index_ktr{k})};
                else
                violin.(viol_met_ktr{k}) = [violin.(viol_met_ktr{k}), ID(i).features.(viol_met_ktr{k})(ID(i).metrics.responder_index_nfkb == 1,viol_met_index_ktr{k})];
                end
            end
           end
        case {'ktr', 'respective'}
           for i = 1:n
            for k = 1:numel(viol_met_ktr)
                if i==1 
                violin.(viol_met_ktr{k}) = {ID(i).features.(viol_met_ktr{k})(ID(i).metrics.responder_index_ktr == 1, viol_met_index_ktr{k})};
                else
                violin.(viol_met_ktr{k}) = [violin.(viol_met_ktr{k}), ID(i).features.(viol_met_ktr{k})(ID(i).metrics.responder_index_ktr == 1,viol_met_index_ktr{k})];
                end
            end
           end
        case 'both'
           for i = 1:n
            for k = 1:numel(viol_met_ktr)
                if i==1 
                violin.(viol_met_ktr{k}) = {ID(i).features.(viol_met_ktr{k})(ID(i).metrics.responder_index_nfkb == 1 & ID(i).metrics.responder_index_ktr == 1,viol_met_index_ktr{k} )};
                else
                violin.(viol_met_ktr{k}) = [violin.(viol_met_ktr{k}), ID(i).features.(viol_met_ktr{k})(ID(i).metrics.responder_index_nfkb == 1 & ID(i).metrics.responder_index_ktr == 1,viol_met_index_ktr{k})];
                end
            end
           end    
    end    

    violin_spacing = 1:n;
    for k = 1:numel(viol_met_ktr)
     %   violin_mack(violin.(viol_met_ktr{k}), violin_spacing, 'Area', 0.05, 'YLim', [-5 15])
        axes.ax(k) = nexttile;
        violin_mack_update(violin.(viol_met_ktr{k}),violin_spacing,'Axes', axes.ax(k), 'Area', 0.15,'XSpace', 0.4, 'BinScale', 1,'Smoothing', 'on', 'Connect', 'on', 'MarkerSize', 7, 'ShowBins', 'off');
    %   violin_kernel(violin.(viol_met_ktr{k}),violin_spacing,'Axes', axes.ax(k), 'Area', 0.04,'XSpace', 0.1, 'BinScale', 1,'Smoothing', 'on', 'Connect', 'on', 'MarkerSize', 7, 'ShowBins', 'off');
        set(axes.ax(k),'XTick',1:n, 'XTickLabel',p.Results.XLabel,'Box','on','GridLineStyle', ':','LineWidth',1,...  
       'YTickLabelMode', 'auto', 'Color', ones(1,3)*1,'XGrid', 'on','GridColor', 'w', 'XTickLabelRotation', 45);
         title([viol_met_ktr{k},' ',num2str(viol_met_index_ktr{k})], 'Interpreter', 'none')
        ylabel(viol_met_units_ktr{k})
    end

end

%%


for i= 1:n
     for k = 1:numel(viol_met_nfkb)
            ID(i).violin.([viol_met_nfkb{k},num2str(viol_met_index_nfkb{k})])= violin.([viol_met_nfkb{k},num2str(viol_met_index_nfkb{k})]){i};
     end
end
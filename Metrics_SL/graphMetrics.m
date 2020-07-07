function [ID] = graphMetrics(IDs, varargin)

%20200623 graphs metrics for Violin Plots, %responders, %oscillators for NFkB and KTR
%TODO make this function graph any features from computeFeatures, not just
%basic metrics (use graphFeatures)

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
addParameter(p, 'StimulationTimePoint', 13, @isnumeric)

%parameter to access metrics to be graphed in violin plots
addParameter(p, 'FeatureListFile', 'C:\Users\stlue\OneDrive\PostDoc UCLA\1 Post Doc UCLA\Matlab analysis\MACKtrack_SL\Metrics_SL\FeatureList.xlsx') %provide file path for Excel table with list of feature to be computed

expectedFilters = {'none','nfkb', 'ktr', 'both', 'respective'};
addParameter(p, 'FilterResponders','none', @(x) any(validatestring(x,expectedFilters)));%filter out non-responders or not

parse(p,IDs, varargin{:})

ViolMetTableNFkB = readtable(p.Results.FeatureListFile, 'Sheet', 'ViolMetTableNFkB');
viol_met_nfkb = table2cell(ViolMetTableNFkB(:,1)); %list of metrics/features to be plotted
viol_met_index_nfkb = table2cell(ViolMetTableNFkB(:,2)); %index for metrics with multiple columns (eg duration, etc.)
viol_met_units_nfkb = table2cell(ViolMetTableNFkB(:,3));
ViolMetTableKTR= readtable(p.Results.FeatureListFile, 'Sheet', 'ViolMetTableKTR');
viol_met_ktr = table2cell(ViolMetTableKTR(:,1)); %list of metrics/features to be plotted
viol_met_index_ktr = table2cell(ViolMetTableKTR(:,2)); %index for metrics with multiple columns (eg duration, etc.)
viol_met_units_ktr = table2cell(ViolMetTableKTR(:,3));

%viol_met_nfkb = {'off_times_nfkb','max_amplitude_nfkb','max_integral_nfkb','pk1_amp_nfkb', 'pk2_amp_nfkb','peakfreq_nfkb','max_derivative_nfkb', 'min_derivative_nfkb', 'pk1_time_nfkb', 'pk2_time_nfkb'};
%viol_met_ktr = {'off_times_ktr','max_amplitude_ktr','max_integral_ktr','pk1_amp_ktr', 'pk2_amp_ktr','peakfreq_ktr','max_derivative_ktr', 'min_derivative_ktr', 'pk1_time_ktr', 'pk2_time_ktr'};

n = numel(IDs);
ID(n).metrics = [];
ID(n).graph = [];
ID(n).info = [];

figure
tiledlayout(4,round(numel(viol_met_nfkb)/2)) %looks best when using even number of metrics/features

% Metrics for NFkB
%run metrics function of desired exp ID
for i= 1:n
    [ID(i).metrics,~,ID(i).graph,ID(i).info,~] = nfkb_ktr_ratio_metrics(IDs(i), 'MinLifetime',p.Results.MinLifetime,...
                            'OnThreshNFkB',p.Results.OnThreshNFkB,'OnThreshKTR',p.Results.OnThreshKTR,...
                            'MinSize', p.Results.MinSize,'StartThreshNFkB', p.Results.StartThreshNFkB,'StartThreshKTR', p.Results.StartThreshKTR, 'Verbose', ... 
                            p.Results.Verbose, 'TrimFrame', p.Results.TrimFrame, 'StimulationTimePoint', p.Results.StimulationTimePoint);
end
%%
%
switch p.Results.FilterResponders 
    case 'none'
       for i = 1:n
        for k = 1:numel(viol_met_nfkb)
            if i==1 
            violin.([viol_met_nfkb{k},num2str(viol_met_index_nfkb{k})]) = {ID(i).metrics.(viol_met_nfkb{k})(:, viol_met_index_nfkb{k})};
            else
            violin.([viol_met_nfkb{k},num2str(viol_met_index_nfkb{k})]) = [violin.([viol_met_nfkb{k},num2str(viol_met_index_nfkb{k})]), ID(i).metrics.(viol_met_nfkb{k})(:, viol_met_index_nfkb{k})];
            end
        end
       end
    case {'nfkb', 'respective'}
       for i = 1:n
        for k = 1:numel(viol_met_nfkb)
            if i==1 
            violin.([viol_met_nfkb{k},num2str(viol_met_index_nfkb{k})]) = {ID(i).metrics.(viol_met_nfkb{k})(ID(i).metrics.responder_index_nfkb == 1, viol_met_index_nfkb{k})};
            else
            violin.([viol_met_nfkb{k},num2str(viol_met_index_nfkb{k})]) = [violin.([viol_met_nfkb{k},num2str(viol_met_index_nfkb{k})]), ID(i).metrics.(viol_met_nfkb{k})(ID(i).metrics.responder_index_nfkb == 1,viol_met_index_nfkb{k})];
            end
        end
       end
    case 'ktr'
       for i = 1:n
        for k = 1:numel(viol_met_nfkb)
            if i==1 
            violin.([viol_met_nfkb{k},num2str(viol_met_index_nfkb{k})]) = {ID(i).metrics.(viol_met_nfkb{k})(ID(i).metrics.responder_index_ktr == 1, viol_met_index_nfkb{k})};
            else
            violin.([viol_met_nfkb{k},num2str(viol_met_index_nfkb{k})]) = [violin.([viol_met_nfkb{k},num2str(viol_met_index_nfkb{k})]), ID(i).metrics.(viol_met_nfkb{k})(ID(i).metrics.responder_index_ktr == 1, viol_met_index_nfkb{k})];
            end
        end
       end
    case 'both'
       for i = 1:n
        for k = 1:numel(viol_met_nfkb)
            if i==1 
            violin.([viol_met_nfkb{k},num2str(viol_met_index_nfkb{k})]) = {ID(i).metrics.(viol_met_nfkb{k})(ID(i).metrics.responder_index_nfkb == 1 & ID(i).metrics.responder_index_ktr == 1, viol_met_index_nfkb{k})};
            else
            violin.([viol_met_nfkb{k},num2str(viol_met_index_nfkb{k})]) = [violin.([viol_met_nfkb{k},num2str(viol_met_index_nfkb{k})]), ID(i).metrics.(viol_met_nfkb{k})(ID(i).metrics.responder_index_nfkb == 1 & ID(i).metrics.responder_index_ktr == 1, viol_met_index_nfkb{k})];
            end
        end
       end   
end

violin_spacing = 1:n;
for k = 1:numel(viol_met_nfkb)
 %   violin_mack(violin.(viol_met_nfkb{k}), violin_spacing, 'Area', 0.05, 'YLim', [-5 15])
    axes.ax(k) = nexttile;
    violin_mack(violin.([viol_met_nfkb{k},num2str(viol_met_index_nfkb{k})]),violin_spacing,'Axes', axes.ax(k), 'Area', 0.04,'XSpace', 0.1, 'BinScale', 1,'Smoothing', 'on', 'Connect', 'on', 'MarkerSize', 7, 'ShowBins', 'off');
%    violin_mack(violin.(viol_met_nfkb{k}),violin_spacing,'Axes', axes.ax(k), 'Area', 0.04,'XSpace', 0.1, 'BinScale', 1,'Smoothing', 'on', 'Connect', 'on', 'MarkerSize', 7, 'ShowBins', 'off');
%    title([viol_met_nfkb{k}], 'Interpreter', 'none')
    title([viol_met_nfkb{k},' ',num2str(viol_met_index_nfkb{k})], 'Interpreter', 'none')
    ylabel(viol_met_units_nfkb{k})
end

%% KTR metrics

%figure
%tiledlayout(2,round(numel(viol_met)/2))

switch p.Results.FilterResponders 
    case 'none'
       for i = 1:n
        for k = 1:numel(viol_met_ktr)
            if i==1 
            violin.([viol_met_ktr{k},num2str(viol_met_index_ktr{k})]) = {ID(i).metrics.(viol_met_ktr{k})(:, viol_met_index_ktr{k})};
            else
            violin.([viol_met_ktr{k},num2str(viol_met_index_ktr{k})]) = [violin.([viol_met_ktr{k},num2str(viol_met_index_ktr{k})]), ID(i).metrics.(viol_met_ktr{k})(:, viol_met_index_ktr{k})];
            end
        end
       end
    case 'nfkb'
       for i = 1:n
        for k = 1:numel(viol_met_ktr)
            if i==1 
            violin.([viol_met_ktr{k},num2str(viol_met_index_ktr{k})]) = {ID(i).metrics.(viol_met_ktr{k})(ID(i).metrics.responder_index_nfkb == 1,viol_met_index_ktr{k})};
            else
            violin.([viol_met_ktr{k},num2str(viol_met_index_ktr{k})]) = [violin.([viol_met_ktr{k},num2str(viol_met_index_ktr{k})]), ID(i).metrics.(viol_met_ktr{k})(ID(i).metrics.responder_index_nfkb == 1,viol_met_index_ktr{k})];
            end
        end
       end
    case {'ktr', 'respective'}
       for i = 1:n
        for k = 1:numel(viol_met_ktr)
            if i==1 
            violin.([viol_met_ktr{k},num2str(viol_met_index_ktr{k})]) = {ID(i).metrics.(viol_met_ktr{k})(ID(i).metrics.responder_index_ktr == 1, viol_met_index_ktr{k})};
            else
            violin.([viol_met_ktr{k},num2str(viol_met_index_ktr{k})]) = [violin.([viol_met_ktr{k},num2str(viol_met_index_ktr{k})]), ID(i).metrics.(viol_met_ktr{k})(ID(i).metrics.responder_index_ktr == 1,viol_met_index_ktr{k})];
            end
        end
       end
    case 'both'
       for i = 1:n
        for k = 1:numel(viol_met_ktr)
            if i==1 
            violin.([viol_met_ktr{k},num2str(viol_met_index_ktr{k})]) = {ID(i).metrics.(viol_met_ktr{k})(ID(i).metrics.responder_index_nfkb == 1 & ID(i).metrics.responder_index_ktr == 1,viol_met_index_ktr{k} )};
            else
            violin.([viol_met_ktr{k},num2str(viol_met_index_ktr{k})]) = [violin.([viol_met_ktr{k},num2str(viol_met_index_ktr{k})]), ID(i).metrics.(viol_met_ktr{k})(ID(i).metrics.responder_index_nfkb == 1 & ID(i).metrics.responder_index_ktr == 1,viol_met_index_ktr{k})];
            end
        end
       end    
end    

violin_spacing = 1:n;
for k = 1:numel(viol_met_ktr)
 %   violin_mack(violin.(viol_met_ktr{k}), violin_spacing, 'Area', 0.05, 'YLim', [-5 15])
    axes.ax(k) = nexttile;
    violin_mack(violin.([viol_met_ktr{k},num2str(viol_met_index_ktr{k})]),violin_spacing,'Axes', axes.ax(k), 'Area', 0.04,'XSpace', 0.1, 'BinScale', 1,'Smoothing', 'on', 'Connect', 'on', 'MarkerSize', 7, 'ShowBins', 'off');
%    violin_mack(violin.(viol_met_ktr{k}),violin_spacing,'Axes', axes.ax(k), 'Area', 0.04,'XSpace', 0.1, 'BinScale', 1,'Smoothing', 'on', 'Connect', 'on', 'MarkerSize', 7, 'ShowBins', 'off');
    title([viol_met_ktr{k},' ',num2str(viol_met_index_ktr{k})], 'Interpreter', 'none')
    ylabel(viol_met_units_ktr{k})
end

%}

%% Percent oscillators plot

osc_cat(n).nfkb = [];
osc_cat(n).ktr = [];
for i = 1:n
osc_cat(i).nfkb    =  get_osc_cats(ID(i).metrics.peakfreq_nfkb,ID(i).metrics.off_times_nfkb,'cutoff_fq', 0.42);
%todo find proper frequency threshold for KTR
osc_cat(i).ktr    =  get_osc_cats(ID(i).metrics.peakfreq_ktr,ID(i).metrics.off_times_ktr,'cutoff_fq', 0.42);            
end

osc_perc_nfkb = nan(1,n);
osc_perc_ktr = nan(1,n);

%calculate %oscillators of total 
switch p.Results.FilterResponders 
    case 'none'
        for i = 1:n
        osc_perc_nfkb(i) = numel(osc_cat(i).nfkb(osc_cat(i).nfkb=='osc'))/numel(osc_cat(i).nfkb);
        osc_perc_ktr(i) = numel(osc_cat(i).ktr(osc_cat(i).ktr=='osc'))/numel(osc_cat(i).ktr);
        end
    case 'both'
        for i = 1:n
        osc_perc_nfkb(i) = numel(osc_cat(i).nfkb(osc_cat(i).nfkb=='osc' & ID(i).metrics.responder_index_nfkb == 1 & ID(i).metrics.responder_index_ktr == 1 ))/numel(osc_cat(i).nfkb((ID(i).metrics.responder_index_nfkb == 1 & ID(i).metrics.responder_index_ktr == 1 )));
        osc_perc_ktr(i) = numel(osc_cat(i).ktr(osc_cat(i).ktr=='osc' & ID(i).metrics.responder_index_nfkb == 1 & ID(i).metrics.responder_index_ktr == 1 ))/numel(osc_cat(i).ktr((ID(i).metrics.responder_index_nfkb == 1 & ID(i).metrics.responder_index_ktr == 1 )));
        end
    case 'ktr'
        for i = 1:n
        osc_perc_nfkb(i) = numel(osc_cat(i).nfkb(osc_cat(i).nfkb=='osc' & ID(i).metrics.responder_index_ktr == 1 ))/numel(osc_cat(i).nfkb( ID(i).metrics.responder_index_ktr == 1 ));
        osc_perc_ktr(i) = numel(osc_cat(i).ktr(osc_cat(i).ktr=='osc' & ID(i).metrics.responder_index_ktr == 1 ))/numel(osc_cat(i).ktr(ID(i).metrics.responder_index_ktr == 1 ));
        end
    case 'nfkb'
        for i = 1:n
        osc_perc_nfkb(i) = numel(osc_cat(i).nfkb(osc_cat(i).nfkb=='osc' & ID(i).metrics.responder_index_nfkb == 1 ))/numel(osc_cat(i).nfkb( ID(i).metrics.responder_index_nfkb == 1 ));
        osc_perc_ktr(i) = numel(osc_cat(i).ktr(osc_cat(i).ktr=='osc' & ID(i).metrics.responder_index_nfkb == 1 ))/numel(osc_cat(i).ktr(ID(i).metrics.responder_index_nfkb == 1 ));
        end
    case 'respective'
        for i = 1:n
        osc_perc_nfkb(i) = numel(osc_cat(i).nfkb(osc_cat(i).nfkb=='osc' & ID(i).metrics.responder_index_nfkb == 1 ))/numel(osc_cat(i).nfkb( ID(i).metrics.responder_index_nfkb == 1 ));
        osc_perc_ktr(i) = numel(osc_cat(i).ktr(osc_cat(i).ktr=='osc' & ID(i).metrics.responder_index_ktr == 1 ))/numel(osc_cat(i).ktr(ID(i).metrics.responder_index_ktr == 1 ));
        end
end       
data_for_osc_bar = [osc_perc_nfkb;osc_perc_ktr];

figure
b = bar(data_for_osc_bar', 'FaceColor', 'flat');
colors = setcolors;
for i = 1:n
    b(1).CData(i,:) = colors.doses{i};
    b(2).CData(i,:) = colors.dosesDarker{i};
end
set(gca, 'XTick', 1:n,'XTickLabels', {IDs});
title('Oscillators [%]')
legend({'NFkB', 'KTR'}, 'Location', 'northwest')
ylabel('Oscillating cells [%]')
xlabel('Experiment ID')



%% Percent responder plot

data_for_resp_nfkb = nan(1,n);
data_for_resp_ktr = nan(1,n);
for i = 1:n
    data_for_resp_nfkb(i) = ID(i).metrics.responders_fraction_nfkb;
    data_for_resp_ktr(i) = ID(i).metrics.responders_fraction_ktr;
end
data_for_resp_bar = [data_for_resp_nfkb;data_for_resp_ktr];
figure
b = bar(data_for_resp_bar', 'FaceColor', 'flat');
colors = setcolors;
for i = 1:n
    b(1).CData(i,:) = colors.doses{i};
    b(2).CData(i,:) = colors.dosesDarker{i};
end

%todo add dose labels as parameter or get from folder name
set(gca, 'XTick', 1:n,'XTickLabels', {IDs});
%set(b, 'DisplayName', {'NFkB', 'KTR'})
title('Responding cells [%]')
legend({'NFkB', 'KTR'}, 'Location', 'northwest')
ylabel('Responder [% of total]')
xlabel('Experiment ID')

%% NFkB vs KTR responders plot
%todo change percent responder plot above to use part of this
responders(n).nfkb = [];
responders(n).ktr = [];
responders(n).dual = [];
responders(n).non = [];
for i = 1:n
    responders(i).nfkb  = numel(ID(i).metrics.responder_index_nfkb(ID(i).metrics.responder_index_nfkb == 1 & ID(i).metrics.responder_index_ktr == 0))/numel(ID(i).metrics.responder_index_nfkb);
    responders(i).ktr   = numel(ID(i).metrics.responder_index_nfkb(ID(i).metrics.responder_index_ktr == 1 & ID(i).metrics.responder_index_nfkb == 0))/numel(ID(i).metrics.responder_index_nfkb);
    responders(i).dual  = numel(ID(i).metrics.responder_index_nfkb(ID(i).metrics.responder_index_ktr == 1 & ID(i).metrics.responder_index_nfkb == 1))/numel(ID(i).metrics.responder_index_nfkb);
    responders(i).non   = numel(ID(i).metrics.responder_index_nfkb(ID(i).metrics.responder_index_ktr == 0 & ID(i).metrics.responder_index_nfkb == 0))/numel(ID(i).metrics.responder_index_nfkb);
end
figure

p1 = plot(1:n,[responders.nfkb], '-o');
p1.Color = [0.9290 0.6940 0.1250];
p1.MarkerEdgeColor = [0.9290 0.6940 0.1250];
p1.MarkerFaceColor = [0.9290 0.6940 0.1250];
hold on
p2 = plot(1:n,[responders.ktr],'-o');
p2.Color = [0, 0.4470, 0.7410];
p2.MarkerEdgeColor = [0, 0.4470, 0.7410];
p2.MarkerFaceColor = [0, 0.4470, 0.7410];
p3 = plot(1:n,[responders.non],'-o');
p3.Color = [113 115 118]/255;
p3.MarkerEdgeColor = [113 115 118]/255;
p3.MarkerFaceColor = [113 115 118]/255;
p4 = plot(1:n,[responders.dual],'-o');
p4.Color = [0.4660 0.6740 0.1880];
p4.MarkerEdgeColor = [0.4660 0.6740 0.1880];
p4.MarkerFaceColor = [0.4660 0.6740 0.1880];

xlim([0.8 (n+0.2)]);
ylim([0 1])
set(gca, 'XTick', 1:n,'XTickLabels', {IDs});
title('Responder Fractions [%]')
legend({'NFkB only', 'KTR only', 'Non-Responders', 'Dual Responders'}, 'Location', 'northeastoutside')
ylabel('Responder Fraction [%]')
xlabel('Experiment ID')
hold off

figure
data_for_resp_cat_bar=[responders(1:n).non;responders(1:n).dual; responders(1:n).nfkb; responders(1:n).ktr]; 
bs = bar(data_for_resp_cat_bar', 'stacked', 'FaceColor', 'flat');
bs(1).FaceColor = [113 115 118]/255;
bs(2).FaceColor =[0.4660 0.6740 0.1880];
bs(3).FaceColor = [0.9290 0.6940 0.1250];
bs(4).FaceColor = [0, 0.4470, 0.7410];
legend({'Non-Responders', 'Dual Responders','NFkB only', 'KTR only'}, 'Location', 'northeastoutside')
title('Responder Fractions [%]')
ylabel('Responder Fraction [%]')
xlabel('Experiment ID')
set(gca, 'XTick', 1:n,'XTickLabels', {IDs});

%}

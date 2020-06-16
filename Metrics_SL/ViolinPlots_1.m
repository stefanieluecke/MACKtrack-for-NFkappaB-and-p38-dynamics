function [ID] = ViolinPlots_1(firstID, lastID)

%20200222 Violin plots for metrics

%run metrics function of desired exp ID

%firstID = 111;
%lastID = 116;

n = lastID - firstID +1;
ID(n).metrics = [];
ID(n).graph = [];
ID(n).info = [];

% Metrics for NFkB
%violin.peakfreq_nfkb = [];
viol_met = {'off_times_nfkb','max_amplitude_nfkb','max_integral_nfkb','pk1_amp_nfkb', 'pk2_amp_nfkb','peakfreq_nfkb','max_derivative_nfkb', 'min_derivative_nfkb', 'pk1_time_nfkb', 'pk2_time_nfkb'};

figure
tiledlayout(4,round(numel(viol_met)/2))
jj = 1;
for i= firstID:lastID
    [ID(jj).metrics,~,ID(jj).graph,ID(jj).info,~] = nfkb_ktr_ratio_metrics(i, "Verbose", 'off', 'MinLifeTime', 100);
jj = jj+1;
end
for j = 1:n
    for k = 1:numel(viol_met)
        if j==1 
        violin.(viol_met{k}) = {ID(j).metrics.(viol_met{k})};
        else
        violin.(viol_met{k}) = [violin.(viol_met{k}), ID(j).metrics.(viol_met{k})];
        end
      end
end

violin_spacing = [1:n];
for k = 1:numel(viol_met)
 %   violin_mack(violin.(viol_met{k}), violin_spacing, 'Area', 0.05, 'YLim', [-5 15])
axes.ax(k) = nexttile;
 violin_mack(violin.(viol_met{k}),violin_spacing,'Axes', axes.ax(k), 'Area', 0.04,'XSpace', 0.1, 'BinScale', 1,'Smoothing', 'on', 'Connect', 'on', 'MarkerSize', 7, 'ShowBins', 'off')
title(viol_met{k}, 'Interpreter', 'none')
end

%% KTR metrics
viol_met = {'off_times_ktr','max_amplitude_ktr','max_integral_ktr','pk1_amp_ktr', 'pk2_amp_ktr','peakfreq_ktr','max_derivative_ktr', 'min_derivative_ktr', 'pk1_time_ktr', 'pk2_time_ktr'};

%figure
%tiledlayout(2,round(numel(viol_met)/2))

for j = 1:n
     for k = 1:numel(viol_met)
        if j==1 
        violin.(viol_met{k}) = {ID(j).metrics.(viol_met{k})};
        else
        violin.(viol_met{k}) = [violin.(viol_met{k}), ID(j).metrics.(viol_met{k})];
        end
      end
 end
violin_spacing = [1:n];
for k = 1:numel(viol_met)
 %   violin_mack(violin.(viol_met{k}), violin_spacing, 'Area', 0.05, 'YLim', [-5 15])
    axes.ax(k) = nexttile;
    violin_mack(violin.(viol_met{k}),violin_spacing,'Axes', axes.ax(k), 'Area', 0.04,'XSpace', 0.1, 'BinScale', 1,'Smoothing', 'on', 'Connect', 'on', 'MarkerSize', 7, 'ShowBins', 'off')
    title(viol_met{k}, 'Interpreter', 'none')
end


%% Percent responder plot
lastID= 116;
firstID= 111;
n = lastID - firstID +1;

data_for_bar_nfkb = nan(1,n);
data_for_bar_ktr = nan(1,n);
for j = 1:n
    data_for_bar_nfkb(j) = ID(j).metrics.responders_fraction_nfkb;
    data_for_bar_ktr(j) = ID(j).metrics.responders_fraction_ktr;
end
data_for_bar = [data_for_bar_nfkb;data_for_bar_ktr];
figure
b = bar(data_for_bar', 'FaceColor', 'flat');
%b = bar(data_for_bar');
colors = setcolors;
for j = 1:n
    b(1).CData(j,:) = colors.doses{j};
    b(2).CData(j,:) = colors.dosesDarker{j};
end
%b(1).DisplayName = {'0', '0.01', '0.1', '1', '10', '100 ng/ml'};
set(gca, 'XTick', 1:n,'XTickLabels', {'0', '0.01', '0.1', '1', '10', '100 ng/ml'});
%set(b, 'DisplayName', {'NFkB', 'KTR'})
title('Responding cells [%]')
legend({'NFkB', 'KTR'}, 'Location', 'northwest')

function [comparisons] = dist_by_JSD(vects)
%INPUT vect: single feature information for all replicates {[rep1_feat_info] ... [repn_feat_info]}
% computes pairwise distance between distributions from hists using
% Jensen–Shannon distance

hists = cell(size(vects));
all = cell2mat(vects(:));
maxval = max(all);  minval = min(all);
%calculate bin width with Freedman Diaconis rule
bin_width = 2*iqr(all)*((numel(all)/length(vects))^(-1/3));
bins = [minval:bin_width:maxval];

for j = 1:length(vects)
    feature_data = vects{j}(~isnan(vects{j}));
%    hists{j} = histcounts(feature_data, bins, 'Normalization', 'pdf');
    hists{j} = histcounts(feature_data, bins, 'Normalization', 'probability');
end

comparisons = nchoosek([1:length(hists)], 2);
comparisons = [comparisons nan(size(comparisons, 1), 1)];
for k = 1:size(comparisons, 1)
comparison_id = comparisons(k, :);
hist1 = hists{comparison_id(1)};
hist2 = hists{comparison_id(2)};
M = 0.5*(hist1+hist2);
JSD = 0;
for i=1:length(hist1)
    if hist1(i) == 0
        add1 = 0;
    else
         add1 = (hist1(i)*bin_width)*(log2((hist1(i)*bin_width)/(M(i)*bin_width))); %% calculate JS divergence
    end
    if hist2(i) == 0
        add2 = 0;
    else
        add2 = (hist2(i)*bin_width)*(log2((hist2(i)*bin_width)/(M(i)*bin_width)));
    end
    JSD = JSD + 0.5*add1 +0.5*add2;
end
comparisons(k, 3) = JSD^0.5;  %take square root of divergence to get distance--true metric
end
end


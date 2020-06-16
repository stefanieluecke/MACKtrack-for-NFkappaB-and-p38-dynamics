function [fold_change, max_fold_change] = get_fold_change(time_series_no_base_ded, baseline, StimulationTimePoint)

%calculates max fold change in time_series matrix (non-baseline deducted)

%fold_change = zeros(size(time_series,1),1);
%max_amp=zeros(size(time_series,1),1);
%min_amp=zeros(size(time_series,1),1);

%max_amp = nanmax(time_series_no_base_ded(:,StimulationTimePoint:end),[],2);
%baseline = nanmean(time_series_no_base_ded(:,1:StimulationTimePoint),2);
fold_change = time_series_no_base_ded(:,StimulationTimePoint:end)./baseline;
max_fold_change = nanmax(fold_change,[],2);

end
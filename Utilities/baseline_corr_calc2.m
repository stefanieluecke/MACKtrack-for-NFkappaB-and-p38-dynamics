%Generating NFkB and KTR Baseline Correction factor

%Provides correction factor for NFkB measurements to correct for drop in
%fluorescence after addition of stimulus media and for decrease in
%fluorescence due to changes in cell morphology over time

% Provides correction factor for KTR measurements to correct for non-stimulus dependent increase
% in signal over time

%Startup 
OneDrivePath = getenv('OneDrive');

% ID# for all mock wells from all data sets to be used in analysis
% Need data up to 8 h ps, thus include dataset which only has 107 TPs almost 9 h total) for MinLifeTime

% Try to generate TP157 long one including all available data for all
% available TPs (latest TPs will be based on less data)
IDs = [111,141,183,117,147,165,129,171,189,135,159,195];
%Redo after picking 2 reps per stimulus to use for final analysis
IDs = [111,141,117,165,129,189,135,159];

[ID] = see_NFkB_KTR_ratio_multiple_HMs(IDs, 'KTRBaselineAdjustment','off','NFkBBaselineAdjustment','off', 'MinLifeTime', 107, 'TrimFrame', 157);


n = 8;
mock_medians_nfkb = [];
mock_medians_ktr = [];
for i = 1:n
    ID(i).metrics.time_series_nfkb(:,end+1:157) = nan;
    ID(i).metrics.time_series_ktr(:,end+1:157) = nan;
    mock_medians_nfkb = [mock_medians_nfkb; nanmedian(ID(i).metrics.time_series_nfkb, 1)];
    mock_medians_ktr = [mock_medians_ktr; nanmedian(ID(i).metrics.time_series_ktr, 1)];
end

NFkBBaselineCorrFact = nanmedian(mock_medians_nfkb,1);
KTRBaselineCorrFact = nanmedian(mock_medians_ktr,1);

save([OneDrivePath, '\PostDoc UCLA\1 Post Doc UCLA\Matlab analysis\MACKtrack_SL\NFkBBaselineAdjustment.mat'], 'NFkBBaselineCorrFact');

save([OneDrivePath, '\PostDoc UCLA\1 Post Doc UCLA\Matlab analysis\MACKtrack_SL\KTRBaselineAdjustment.mat'], 'KTRBaselineCorrFact');

%% Testing the result
%{
[ID3] = see_NFkB_KTR_ratio_multiple_HMs(IDs, 'KTRBaselineAdjustment','off','NFkBBaselineAdjustment','off', 'MinLifeTime', 107, 'TrimFrame', 107, 'SortMetric', 'max_amplitude_ktr');

[ID4] = see_NFkB_KTR_ratio_multiple_HMs(111:116, 'KTRBaselineAdjustment','on','NFkBBaselineAdjustment','on', 'MinLifeTime', 107, 'TrimFrame', 107, 'SortMetric', 'max_amplitude_ktr');
[ID5] = see_NFkB_KTR_ratio_multiple_HMs(111:116, 'KTRBaselineAdjustment','off','NFkBBaselineAdjustment','off', 'MinLifeTime', 107, 'TrimFrame', 107, 'SortMetric', 'max_amplitude_ktr');

[ID2] = see_NFkB_KTR_ratio_multiple_means(IDs, 'KTRBaselineAdjustment','on','NFkBBaselineAdjustment','on', 'MinLifeTime', 107, 'TrimFrame', 107, 'GraphLimitsKTR', [-0.1 0.35],'GraphLimitsNFkB', [-2 7]);
[ID2] = see_NFkB_KTR_ratio_multiple_means(IDs(1:6), 'KTRBaselineAdjustment','on','NFkBBaselineAdjustment','on', 'MinLifeTime', 107, 'TrimFrame', 107, 'GraphLimitsKTR', [-0.1 0.35],'GraphLimitsNFkB', [-2 7]);
[ID2] = see_NFkB_KTR_ratio_multiple_means(IDs(7:12), 'KTRBaselineAdjustment','on','NFkBBaselineAdjustment','on', 'MinLifeTime', 107, 'TrimFrame', 107);

[ID3] = see_NFkB_KTR_ratio_multiple_means(IDs, 'KTRBaselineAdjustment','off','NFkBBaselineAdjustment','off', 'MinLifeTime', 107, 'TrimFrame', 107, 'GraphLimitsKTR', [-0.1 0.35],'GraphLimitsNFkB', [-2 7]);
[ID3] = see_NFkB_KTR_ratio_multiple_means(IDs(1:6), 'KTRBaselineAdjustment','off','NFkBBaselineAdjustment','off', 'MinLifeTime', 107, 'TrimFrame', 107, 'GraphLimitsKTR', [-0.1 0.35],'GraphLimitsNFkB', [-2 7]);
[ID3] = see_NFkB_KTR_ratio_multiple_means(IDs(7:12), 'KTRBaselineAdjustment','off','NFkBBaselineAdjustment','off', 'MinLifeTime', 107, 'TrimFrame', 107);


[ID4] = see_NFkB_KTR_ratio_multiple_means(111:116, 'KTRBaselineAdjustment','on','NFkBBaselineAdjustment','on', 'MinLifeTime', 107, 'TrimFrame', 107);
[ID5] = see_NFkB_KTR_ratio_multiple_means(111:116, 'KTRBaselineAdjustment','off','NFkBBaselineAdjustment','off', 'MinLifeTime', 107, 'TrimFrame', 107);
%}

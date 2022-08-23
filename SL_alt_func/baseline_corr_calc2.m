%Generating NFkB and KTR Baseline Correction factor

%Provides correction factor for NFkB measurements to correct for drop in
%fluorescence after addition of stimulus media and for decrease in
%fluorescence due to changes in cell morphology over time

% Provides correction factor for KTR measurements to correct for non-stimulus dependent increase
% in signal over time

load('C:\Users\stlue\OneDrive\PostDoc UCLA\1 Post Doc UCLA\Matlab analysis\testing\20200706 Testing additional baseline corrections\ID_mock_collection2.mat')

IDs = [];
[ID] = see_NFkB_KTR_ratio_multiple_HMs(IDs, ,'KTRBaselineAdjustment','off','NFkBBaselineAdjustment','off');


n = 12;
mock_medians_nfkb = [];
mock_medians_ktr = [];
for i = 1:n
    mock_medians_nfkb = [mock_medians_nfkb; nanmedian(ID(i).metrics.time_series_nfkb, 1)];
    mock_medians_ktr = [mock_medians_ktr; nanmedian(ID(i).metrics.time_series_ktr, 1)];
end

NFkBBaselineCorrFact = 0 - nanmedian(mock_medians_nfkb,1);
KTRBaselineCorrFact = 0 - nanmedian(mock_medians_ktr,1);

save('NFkBBaselineAdjustment.mat', 'NFkBBaselineCorrFact');

save('KTRBaselineAdjustment.mat', 'KTRBaselineCorrFact');

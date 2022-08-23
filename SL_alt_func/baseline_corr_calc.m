%Generating NFkB Baseline Correction factor

%Provides correction factor for NFkB measurements to correct for drop in
%fluorescence after addition of stimulus media and for decrease in
%fluorescence due to changes in cell morphology over time

load('C:\Users\stlue\OneDrive\PostDoc UCLA\1 Post Doc UCLA\Matlab analysis\testing\20200706 Testing additional baseline corrections\ID_mock_collection2.mat')



n = 7;
mock_medians = [];
for i = 1:n
    mock_medians = [mock_medians; nanmedian(ID(i).metrics.time_series_nfkb, 1)];
end

NFkBBaselineCorrFact = 0 - nanmedian(mock_medians,1);
save('NFkBBaselineAdjustment.mat', 'NFkBBaselineCorrFact');

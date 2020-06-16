function [HMW,Prominence, Amp, Locs] = halfMaxWidth(time_series, PeakHeight, FPH, StimulationTimePoint)
%--------------------------------------------------------------------------
%[HMW,Prominence, Amp, Locs] = halfMaxWidth(time_series, PeakHeight)
%--------------------------------------------------------------------------
%todo adjust to new StimulationTimePoint!
% if isrow(PeakHeight)
%    PeakHeight = transpose( PeakHeight);
% end 7/10/2019

eqDims= size(time_series(:,StimulationTimePoint:end)) == size(PeakHeight); 
if ~eqDims(1) && eqDims(2)
    %check that dimensions of time_series and PeakHeights are aligned
     PeakHeight = PeakHeight.' ;
end
nPeaks=size(PeakHeight,2);
HMW = nan(size(PeakHeight,1), nPeaks);

Amp=HMW; Locs = HMW; Prominence = HMW; 
lastFrame =size(time_series(:,StimulationTimePoint:end),2);
t=linspace(0,(lastFrame-1)/FPH, lastFrame);

%todo check these parameters, etc
%todo this gives error for low stimulations
    for i =1:size(time_series(:,StimulationTimePoint:end),1)
        minPeakHgt=min(PeakHeight(i,:)*0.95);
        if isnan(minPeakHgt); continue; end
    %20200608 test because low stim are giving problem here     
          %  [Amp(i,:), Locs(i,:), HMW(i,:), Prominence(i,:)] =findpeaks(time_series(i,StimulationTimePoint:end), t,'NPeaks',nPeaks,...
           % 'MinPeakHeight', minPeakHgt,'WidthReference', 'halfheight', 'Annotate', 'extents', 'MinPeakDistance', 6/FPH);
   %TODO check the effect of relacing WidthReference halfheight with default halfProminence (halfheight gave problems with peaks below 0)
 %TODO also check effect of removing MinPeakHeight requirement, because
 %this gave problem iwth mock conditons (due to neg number due to new baseline and *0.95PeakHeight determination of minPeakHeight
            [Amp(i,:), Locs(i,:), HMW(i,:), Prominence(i,:)] =findpeaks(time_series(i,StimulationTimePoint:end), t,'NPeaks',nPeaks,...
            'Annotate', 'extents', 'MinPeakDistance', 6/FPH);

    end
end
function  time_to_half_max_integral=get_time_to_half_max_integral(integrals, FramesPerHour, StimulationTimePoint)
%note: in the callings of this function in computeFeatures and
%get_pos_integrals, 1 is used as StimulationTimePoint, because integrals
%were cut to only be calculated after StimulationTimePoint

halfMaxIntegral = nanmax(integrals(:,StimulationTimePoint:end),[],2)/2;
 
 distances = abs(integrals(:,StimulationTimePoint:end)- halfMaxIntegral);
 [~, idx] = nanmin(distances,[],2);
 idx(idx==1) = NaN;
 time_to_half_max_integral = (idx-1)/FramesPerHour;

end
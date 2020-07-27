function [pos_integral_features_ktr] = get_pos_integrals_ktr(integrals, FramesPerHour, endFrame)

%integrals now already are calculated using only values after StStimulationTimePoint 
%integrals = integrals(:, StimulationTimePoint:end);

integrals_pos = nan(size(integrals));
integrals_pos(:,1) = integrals(:,1);
for i = 1:size(integrals_pos,1)
    for j = 2:size(integrals_pos,2)
        if integrals(i,j)>= integrals(i,j-1)
         integrals_pos(i,j) =  integrals_pos(i,j-1) + (integrals(i,j) - integrals(i,j-1));
        %elseif integrals(i,j) = integrals(i,j-1)
        elseif integrals(i,j)< integrals(i,j-1)
         integrals_pos(i,j) =  integrals_pos(i,j-1);
        end
    end
end

pos_integral_features_ktr.integrals_pos_ktr =integrals_pos; 

pos_integral_features_ktr.max_pos_integral_ktr = nanmax(integrals_pos,[],2);

%here '1' instead of 'StimulationTimePoint', because size of pos integrals already only includes values from StimulationTimePoint onwards
pos_integral_features_ktr.time2HalfMaxPosIntegral_ktr = get_time_to_half_max_integral(integrals_pos(:,1:endFrame), FramesPerHour, 1);
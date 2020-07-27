function [pos_integral_features_nfkb] = get_pos_integrals_nfkb(integrals, FramesPerHour, endFrame)

%integrals now are already calculated using only StimulationTimePoint onwards 
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

pos_integral_features_nfkb.integrals_pos_nfkb =integrals_pos; 

pos_integral_features_nfkb.max_pos_integral_nfkb = nanmax(integrals_pos,[],2);

%here '1' instead of 'StimulationTimePoint', because size of pos integrals already only includes values from StimulationTimePoint onwards
pos_integral_features_nfkb.time2HalfMaxPosIntegral_nfkb = get_time_to_half_max_integral(integrals_pos(:,1:endFrame), FramesPerHour, 1);
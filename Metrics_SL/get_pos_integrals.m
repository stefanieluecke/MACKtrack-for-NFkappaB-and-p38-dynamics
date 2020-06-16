function [pos_integral_features] = get_pos_integrals(integrals, StimulationTimePoint)

%integrals are shortened to include only values from StimulationTimePoint onwards
integrals = integrals(:, StimulationTimePoint:end);

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

pos_integral_features.integrals_pos

 = get_time_to_half_max_integral(metrics.integrals_nfkb(:,1:endFrame), FramesPerHour, StimulationTimePoint);
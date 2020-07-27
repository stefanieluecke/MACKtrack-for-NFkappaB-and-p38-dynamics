function [max_pk1_speed,max_pk1_speed_time] = get_max_pk1_speed(pk1_time, derivatives, FramesPerHour, StimulationTimePoint)

%todo check whether the smoothing of derivatives is approppriate, Ade uses
%different smoothing function that I can't find in his Github
smoothed_derivatives = smoothdata(derivatives,'lowess'); 
pk1_frame = pk1_time *FramesPerHour + StimulationTimePoint;
max_pk1_speed = zeros(size(pk1_time));
max_pk1_speed_frame = zeros(size(pk1_time));
for i = 1:numel(max_pk1_speed)
    if isnan(pk1_time(i))
        max_pk1_speed(i)   = NaN;
        max_pk1_speed_frame(i) = NaN;
    else
        [max_pk1_speed(i), max_pk1_speed_frame(i)] = nanmax(smoothed_derivatives(i,StimulationTimePoint:pk1_frame(i)),[],2);

    end
end
max_pk1_speed_time = max_pk1_speed_frame./FramesPerHour;
end
function [max_pk1_speed,max_pk1_speed_frame] = get_max_pk1_speed(pk1_time, derivatives, FramesPerHour, StimulationTimePoint)

%todo check whether the smoothing of derivatives is approppriate, Ade uses
%different smoothing function that I can't find in his Github
smoothed_derivatives = smoothrows(derivatives,3); 
pk1_frame = pk1_time *FramesPerHour + StimulationTimePoint;
for 
    if pk1_time == NaN
    max_pk1_speed   = NaN
    
    [max_pk1_speed, max_pk1_speed_frame] = nanmax(smoothed_derivatives(:,StimulationTimePoint:pk1_frame),[],2);

end
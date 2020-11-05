function output= get_peak_stats_ktr(time_series,baseline_stdv, varargin)

%Calculates peak statistics
%--------------------------------------------------------------------------

p=inputParser;
%addRequired(p,'time_series', @isstruct);
%addRequired(p,'pk1_amp', @isstruct);
%addRequired(p,'pk2_amp', @isstruct);

addRequired(p,'time_series');
addRequired(p,'baseline_stdv');

addParameter (p,'min_pks',2,@isnumeric);
addParameter (p,'max_pk_diff',35,@isnumeric);
addParameter(p, 'StimulationTimePoint', 13, @isnumeric);
addParameter(p, 'FramesPerHour', 12, @isnumeric);
parse (p,time_series, baseline_stdv,varargin{:});
%time_series = p.Results.time_series;
%pk1_amp = p.Results.pk1_amp;
%pk2_amp = p.Results.pk2_amp;
min_pks =p.Results.min_pks;
max_pk_diff=p.Results.max_pk_diff;
FramesPerHour = p.Results.FramesPerHour;
StimulationTimePoint = p.Results.StimulationTimePoint;

%% Peak calling and filtering
%todo adjust nfkbpeaks parameters for KTR
%todo adjust these for new starting timepoint!!!
%todo check for any differences between Ade's and Brooks' nfkbpeaks
%[output.peak_times_ktr,output.peak_amps_ktr, output.valley_times_ktr, output.valley_amps_ktr]=nfkbpeaks(time_series, 'BeginFrame',3,'MinHeight',0.75,'MinDist',6);
[output.peak_times_ktr,output.peak_amps_ktr, output.valley_times_ktr, output.valley_amps_ktr]=nfkbpeaks(time_series(:, StimulationTimePoint:end),baseline_stdv, 'BeginFrame',3,'MinHeight',0.03,'MinDist',4, 'SmoothSize',1);

%Smoothing --> incl in nfkbpeaks, not in metrics
%amount of peaks to ask globalpeaks to make
%get valley times, valley amps

%% Calculating peak features
ipt=diff(output.peak_times_ktr,1,2);

ipt(ipt>max_pk_diff)=nan;

tot_pks = sum(~isnan(ipt),2)+1;
output.kept_ktr = true([size(ipt,1),1]);

ipt(tot_pks<min_pks,:) = [];
output.kept_ktr(tot_pks<min_pks)=false;
ipt = ipt.*5;%convert to minutes

%Brooks' analysis excludes ipt > 35, excludes cells with <5peaks

output.mean_ipt_ktr    =nanmean(ipt,2);
output.median_ipt_ktr  =nanmedian(ipt,2);
output.std_ipt_ktr     =nanstd(ipt,[],2);
% output.var_ipt        =nanvar(ipt,[],2);
output.max_ipt_ktr     = nanmax(ipt,[],2);
output.min_ipt_ktr     = nanmin(ipt,[],2);
output.cv_ipt_ktr      = output.std_ipt_ktr./output.mean_ipt_ktr;
output.ipt_ktr         =ipt;
output.num_peaks_ktr   =tot_pks;
%%
output.mean_peak_amp_ktr   = nanmean(output.peak_amps_ktr, 2); 
output.median_peak_amp_ktr = nanmedian(output.peak_amps_ktr, 2); 
% output.var_peak_amp_ktr  = nanvar(output.peak_amps_ktr, [],2); 
output.std_peak_amp_ktr    = nanstd(output.peak_amps_ktr, [],2); 
output.cv_peak_amp_ktr     = output.std_peak_amp_ktr./output.mean_peak_amp_ktr; 

% output.peak2trough =nan(size(output.valley_amps,1) , size(output.valley_amps,2)*2); 
% peak2trough =output.peak_amps(:,1:end-1)-output.valley_amps; 
% output.peak2trough(:,1:2:end) = peak2trough; 
% output.peak2trough(:, 2:2:end) = output.peak_amps(:,2:end)-output.valley_amps;


%todo understand the calculation of these metrics!

output.peak2trough_ktr     = output.peak_amps_ktr(:,1:end-1)-output.valley_amps_ktr;
minVals = zeros(size(time_series,1),1);
for row=1:size(minVals,1)
    pk_frame =output.peak_times_ktr(row, 1);
    if isnan(pk_frame)
        continue;
    end
    minVals(row) = min(time_series(row, StimulationTimePoint:pk_frame+StimulationTimePoint));
end
troughs = [minVals, output.valley_amps_ktr];

output.trough2peak_ktr     = troughs-output.peak_amps_ktr;

% output.peak2
output.max_peak2trough_ktr     =nanmax(output.peak2trough_ktr, [],2);
output.max_trough2peak_ktr     =nanmax(output.trough2peak_ktr, [],2);

output.mean_peak2trough_ktr    =nanmean(output.peak2trough_ktr, 2);
output.mean_trough2peak_ktr    =nanmean(output.trough2peak_ktr,2);

output.median_peak2trough_ktr  =nanmedian(output.peak2trough_ktr,2); 
output.median_trough2peak_ktr  =nanmedian(output.trough2peak_ktr,2);

output.min_peak2trough_ktr     =nanmin(output.peak2trough_ktr,[],2); 
output.min_trough2peak_ktr     =nanmin(output.trough2peak_ktr, [],2);

output.std_trough2peak_ktr    =nanstd(output.trough2peak_ktr,[],2);
output.std_peak2trough_ktr    =nanstd(output.peak2trough_ktr,[],2); 
output.cv_trough2peak_ktr     =output.std_trough2peak_ktr./output.mean_trough2peak_ktr;
output.cv_peak2trough_ktr     =output.std_peak2trough_ktr./output.mean_peak2trough_ktr; 

%{
%todo compare globalpeaks used in this to the one used above
pkAmps = [pk1_amp,pk2_amp];
[pkWidth,pkProminence,~,~] = halfMaxWidth(time_series,pkAmps, FramesPerHour, StimulationTimePoint);
output.pk1_width_ktr =pkWidth(:,1); 
output.pk2_width_ktr =pkWidth(:,2); 
output.pk1_prom_ktr = pkProminence(:,1); 
output.pk2_prom_ktr = pkProminence(:,2); 
%}
end
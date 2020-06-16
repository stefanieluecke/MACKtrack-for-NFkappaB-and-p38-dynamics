function output= get_peak_stats_nfkb(time_series,pk1_amp, pk2_amp, varargin)

%Calculates peak statistics
%--------------------------------------------------------------------------

%todo double check thoroughly that StimulationTimePoint adjustment worked
%todo why does IPT have less rows than number of cells

p=inputParser;
%addRequired(p,'time_series', @isstruct);
%addRequired(p,'pk1_amp', @isstruct);
%addRequired(p,'pk2_amp', @isstruct);

addRequired(p,'time_series');
addRequired(p,'pk1_amp');
addRequired(p,'pk2_amp');

addParameter (p,'min_pks',2,@isnumeric);
addParameter (p,'max_pk_diff',35,@isnumeric);
addParameter(p, 'StimulationTimePoint', 13, @isnumeric);
addParameter(p, 'FramesPerHour', 12, @isnumeric);
parse (p,time_series, pk1_amp, pk2_amp,varargin{:});
%time_series = p.Results.time_series;
%pk1_amp = p.Results.pk1_amp;
%pk2_amp = p.Results.pk2_amp;
min_pks =p.Results.min_pks;
max_pk_diff=p.Results.max_pk_diff;
FramesPerHour = p.Results.FramesPerHour;
StimulationTimePoint = p.Results.StimulationTimePoint;

%todo adjust nfkbpeaks for KTR
%todo adjust these for new starting timepoint
%todo check for any differences between Ade's and Brooks' nfkbpeaks
%todo check fi these parameters work for me
[output.peak_times_nfkb,output.peak_amps_nfkb, output.valley_times_nfkb, output.valley_amps_nfkb]=nfkbpeaks(time_series(:, StimulationTimePoint:end), 'BeginFrame',3,'MinHeight',0.75,'MinDist',6);
ipt=diff(output.peak_times_nfkb,1,2);

%%
ipt(ipt>max_pk_diff)=nan;

tot_pks = sum(~isnan(ipt),2)+1;
output.kept_nfkb = true([size(ipt,1),1]);

ipt(tot_pks<min_pks,:) = [];
output.kept_nfkb(tot_pks<min_pks)=false;
ipt = ipt.*5;%convert to minutes

%Brooks' analysis excludes ipt > 35, excludes cells with <5peaks

output.mean_ipt_nfkb    =nanmean(ipt,2);
output.median_ipt_nfkb  =nanmedian(ipt,2);
output.std_ipt_nfkb     =nanstd(ipt,[],2);
% output.var_ipt        =nanvar(ipt,[],2);
output.max_ipt_nfkb     = nanmax(ipt,[],2);
output.min_ipt_nfkb     = nanmin(ipt,[],2);
output.cv_ipt_nfkb      = output.std_ipt_nfkb./output.mean_ipt_nfkb;
output.ipt_nfkb         =ipt;
output.num_peaks_nfkb   =tot_pks;
%%
output.mean_peak_amp_nfkb   = nanmean(output.peak_amps_nfkb, 2); 
output.median_peak_amp_nfkb = nanmedian(output.peak_amps_nfkb, 2); 
% output.var_peak_amp_nfkb  = nanvar(output.peak_amps_nfkb, [],2); 
output.std_peak_amp_nfkb    = nanstd(output.peak_amps_nfkb, [],2); 
output.cv_peak_amp_nfkb     = output.std_peak_amp_nfkb./output.mean_peak_amp_nfkb; 

% output.peak2trough =nan(size(output.valley_amps,1) , size(output.valley_amps,2)*2); 
% peak2trough =output.peak_amps(:,1:end-1)-output.valley_amps; 
% output.peak2trough(:,1:2:end) = peak2trough; 
% output.peak2trough(:, 2:2:end) = output.peak_amps(:,2:end)-output.valley_amps;

output.peak2trough_nfkb     = output.peak_amps_nfkb(:,1:end-1)-output.valley_amps_nfkb;
minVals = zeros(size(time_series,1),1);
for row=1:size(minVals,1)
    pk_frame =output.peak_times_nfkb(row, 1);
    if isnan(pk_frame)
        continue;
    end
    minVals(row) = min(time_series(row, StimulationTimePoint:pk_frame+StimulationTimePoint));
end
troughs = [minVals, output.valley_amps_nfkb];

output.trough2peak_nfkb     = troughs-output.peak_amps_nfkb;

% output.peak2
output.max_peak2trough_nfkb     =nanmax(output.peak2trough_nfkb, [],2);
output.max_trough2peak_nfkb     =nanmax(output.trough2peak_nfkb, [],2);

output.mean_peak2trough_nfkb    =nanmean(output.peak2trough_nfkb, 2);
output.mean_trough2peak_nfkb    =nanmean(output.trough2peak_nfkb,2);

output.median_peak2trough_nfkb  =nanmedian(output.peak2trough_nfkb,2); 
output.median_trough2peak_nfkb  =nanmedian(output.trough2peak_nfkb,2);

output.min_peak2trough_nfkb     =nanmin(output.peak2trough_nfkb,[],2); 
output.min_trough2peak_nfkb     =nanmin(output.trough2peak_nfkb, [],2);

output.std_trough2peak_nfkb    =nanstd(output.trough2peak_nfkb,[],2);
output.std_peak2trough_nfkb    =nanstd(output.peak2trough_nfkb,[],2); 
output.cv_trough2peak_nfkb     =output.std_trough2peak_nfkb./output.mean_trough2peak_nfkb;
output.cv_peak2trough_nfkb     =output.std_peak2trough_nfkb./output.mean_peak2trough_nfkb; 

%todo compare globalpeaks used in this to the one used above
pkAmps = [pk1_amp,pk2_amp];
[pkWidth,pkProminence,~,~] = halfMaxWidth(time_series,pkAmps, FramesPerHour, StimulationTimePoint);
output.pk1_width_nfkb =pkWidth(:,1); 
output.pk2_width_nfkb =pkWidth(:,2); 
output.pk1_prom_nfkb = pkProminence(:,1); 
output.pk2_prom_nfkb = pkProminence(:,2); 

end




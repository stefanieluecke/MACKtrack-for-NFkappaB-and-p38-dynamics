function [metrics,aux, graph, info, measure] = nfkb_ktr_nuc_metrics(id,varargin)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% metrics = nfkb_ktr_metrics(id)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% NFKBMETRICS uses the see_nfkb_native function to filter and preprocess NFkB trajectories,
% then calculates related metrics regarding activation. Metric types include:
% 
% 1) time series (base NFkB and ktr dynamics, resampled to 12 frames/hr
% 2) integrated activity
% 3) differentiated activity
% 4) calculated metrics: measuring aspects of oscillation, duration, timing ,and amplitude
%
% INPUTS (required):
% id             filename or experiment ID (from Google Spreadsheet specified in "locations.mat")
%
% INPUT PARAMETERS (optional; specify with name-value pairs)
% 'Display'         'on' or 'off' - show graphs (default: process data only; no graphs)
% 'Verbose'          'on' or 'off' - show verbose output
% 'MinLifetime'      final frame used to filter for long-lived cells (default = 100)
% 'TrimFrame'        trim sets to common length (default = 254 timepoints) 
% 'ConvectionShift'  max allowed time shift between scenes (to correct for poor mixing - default is no shift allowed)
% 'MinSize'
%
% OUTPUT: 
% metrics   structure with output fields
% aux       Extra data (e.g. fourier information (FFT, power, frequencies), thresholds used in envelope/duration)
% graph     main structure output from see_nfkb_native
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
%% INPUT PARSING
% Create input parser object, add required params from function input
p = inputParser;
% Required: ID input
valid_id = @(x) assert((isnumeric(x)&&length(x)==1)||exist(x,'file'),...
    'ID input must be spreadsheet ID or full file path');
addRequired(p,'id',valid_id);
% Optional parameters
addParameter (p, 'OnThreshNFkB', 2, @isnumeric);
addParameter (p, 'OnThreshKTR', 2, @isnumeric);
addParameter(p, 'StartThreshNFkB', 10, @isnumeric); 
addParameter(p,'StartThreshKTR',8, @isnumeric);
addParameter(p, 'MinSize', 90, @isnumeric); 
addParameter(p,'MinLifetime',117, @isnumeric);
addParameter(p,'TrimFrame',157, @isnumeric);
addParameter (p, 'GraphLimitsNFkB',[-0.25 8],@isnumeric);
addParameter (p, 'GraphLimitsKTR',[0 500],@isnumeric);
expectedFlags = {'on','off'};
addParameter(p,'Verbose','off', @(x) any(validatestring(x,expectedFlags)))
valid_conv = @(x) assert(isnumeric(x)&&(x>=0)&&(length(x)==1),...
    'Convection correction parameter must be single integer >= 0');
addParameter(p,'ConvectionShift',1, valid_conv);
addParameter(p, 'StartTimePoint', 13, @isnumeric)

parse(p,id, varargin{:})

%% PARAMETETERS for finding off times - chosen using 'scan_off_params.m'
OnThreshNFkB = p.Results.OnThreshNFkB; % Minimum activity required for cell to register as 'on'
OnThreshKTR = p.Results.OnThreshKTR; % Minimum activity required for cell to register as 'on'
window_sz = 14; % ~1 hr windows (on either side of a given timepoint)
thresh = 0.9; % Pct of inactivity allowed in a given window
cutoff_time = 4; % time to look for cell activity before declaring it "off" (hrs)
off_pad = 12; % Signal time added to trajectory in  FFT calculation (keeps transients from being recorded as osc.)

%% INITIALIZATION. Load and process data. Interpolate time series, calculate deriv/integral approximations
%{
%Brooks' version commented out, I'm trying Ade's version (below)
if ~ismember('MinLifetime',p.UsingDefaults)
   [graph, info, measure] = see_nfkb_native(id,'MinLifetime',p.Results.MinLifetime,...
                            'ConvectionShift',p.Results.ConvectionShift);
   graph.var = graph.var(:,1:p.Results.MinLifetime);
   graph.t = graph.t(1:size(graph.var,2));
else
   [graph, info, measure] = see_nfkb_native(id, 'ConvectionShift',p.Results.ConvectionShift);
end
%}

StartThreshNFkB = p.Results.StartThreshNFkB; 
StartThreshKTR = p.Results.StartThreshKTR; 
MinSize = p.Results.MinSize; 
MinLifetime = p.Results.MinLifetime; 
ConvectionShift = p.Results.ConvectionShift; 


%? double check all parameter passing through makes sense, esp baseline!
[graph, info, measure] = filter_nfkb_ktr_nuc(id,'MinLifetime',MinLifetime,...
                            'ConvectionShift',ConvectionShift, 'OnThreshNFkB',OnThreshNFkB,...
                            'OnThreshKTR',OnThreshKTR,'MinSize', MinSize,'StartThreshNFkB',...
                            StartThreshNFkB,'StartThreshKTR', StartThreshKTR, 'Verbose', p.Results.Verbose,...
                            'GraphLimitsNFkB', p.Results.GraphLimitsNFkB, 'GraphLimitsKTR', p.Results.GraphLimitsKTR, 'StartTimePoint', p.Results.StartTimePoint);
   
graph.var_nfkb = graph.var_nfkb(:,1:min(p.Results.TrimFrame, size(graph.var_nfkb,2)));
graph.var_ktr = graph.var_ktr(:,1:min(p.Results.TrimFrame, size(graph.var_ktr,2)));
graph.t = graph.t(1:size(graph.var_nfkb,2));

graph.opt_nfkb = maketicks(graph.t,info.GraphLimitsNFkB,0); %calls function to add tick labels to time frame (?only partially understand what that does)
graph.opt_ktr = maketicks(graph.t,info.GraphLimitsKTR,0); %calls function to add tick labels to time frame (?only partially understand what that does)
graph.opt_nfkb.Name = 'NF\kappaB activation'; 
graph.opt_ktr.Name = 'Kinase activation'; 


if ~ismember ('OnThreshNFkB',p.UsingDefaults)
    OnThreshNFkB = info.OnThreshNFkB;
end
if ~ismember ('OnThreshKTR',p.UsingDefaults)
    OnThreshKTR = info.OnThreshKTR;
end

%% NFkB METRICS
%% BASIC NFkB METRICS: TIME SERIES, DERIVATIVE, INTEGRAL
% 1) basic time series. Interpolate over "normal" interval (12 frames per hr) if required
t = min(graph.t):1/12:max(graph.t);
if length(t)~=length(graph.t)
    metrics.time_series_nfkb = nan(size(graph.var_nfkb,1),length(t));
    for i = 1:size(graph.var_nfkb,1)
        metrics.time_series_nfkb(i,:) = interp1(graph.t,graph.var_nfkb(i,:),t);
    end
else
    metrics.time_series_nfkb = graph.var_nfkb;
end

% 2) integrated activity
metrics.integrals_nfkb = nan(size(metrics.time_series_nfkb));
for i = 1:size(metrics.integrals_nfkb,1)
    metrics.integrals_nfkb(i,:) = cumtrapz(t,metrics.time_series_nfkb(i,:));
end

% 3) differentiated activity - use central finite difference
smoothed = smoothrows(metrics.time_series_nfkb,3);
metrics.derivatives_nfkb = (smoothed(:,3:end) - smoothed(:,1:end-2))/(1/6);


%% TRIM EVERYBODY to a common length (of "good" sets, current minimum is roughly 21 hrs)
try
    metrics.time_series_nfkb = metrics.time_series_nfkb(:,1:p.Results.TrimFrame);
    metrics.integrals_nfkb = metrics.integrals_nfkb(:,1:p.Results.TrimFrame);
    metrics.derivatives_nfkb = metrics.derivatives_nfkb(:,1:(p.Results.TrimFrame-2));
    smoothed = smoothed(:,1:p.Results.TrimFrame);
    t = t(1:p.Results.TrimFrame);

catch me
    disp(['Note: vectors too short to cap @ ',num2str(p.Results.TrimFrame),' frames'])
end
%% MISC NFkB METRICS

% 4) Integrals within one-hour windows (0-1, 1-2, 2-3) and three hour windows (0-3, 1-4, etc) of activity
max_hr = floor(max(t));
metrics.intwin1_nfkb = nan(size(metrics.time_series_nfkb,1),max_hr);
metrics.intwin3_nfkb = nan(size(metrics.time_series_nfkb,1),max_hr-2);
for i = 1:(max_hr)
    win = t>=(i-1) & t<(i);
    metrics.intwin1_nfkb(:,i) = trapz(t(win),metrics.time_series_nfkb(:,win),2);
    if i<= (max_hr-2)
        win = t>=(i-1) & t<(i+2);
        metrics.intwin3_nfkb(:,i) = trapz(t(win),metrics.time_series_nfkb(:,win),2);
    end
end

% MAX/MIN metrics
metrics.max_amplitude_nfkb = nanmax(metrics.time_series_nfkb,[],2);
metrics.max_integral_nfkb = nanmax(metrics.integrals_nfkb,[],2);
metrics.max_derivative_nfkb = nanmax(metrics.derivatives_nfkb,[],2);
metrics.min_derivative_nfkb = nanmin(metrics.derivatives_nfkb,[],2);

% ACTIVITY metrics: compute an off-time for all cells
metrics.off_times_nfkb = zeros(size(smoothed,1),1);
inactive = [repmat(nanmin(smoothed(:,1:7),[],2),1,window_sz*2+1),smoothed(:,:),...
    repmat(nanmedian(smoothed(:,(end-window_sz:end)),2),1,window_sz*2)];
inactive = smoothrows(inactive<(OnThreshNFkB),(window_sz*2));
frontcrop = round(window_sz*2*(1-thresh))+window_sz+1;
inactive = inactive(:,frontcrop:end);
inactive = inactive(:,1:size(smoothed,2));
inactive(isnan(smoothed)) = nan;

% Find the final time each cell was active
for i = 1:length(metrics.off_times_nfkb)
    active_times = find(inactive(i,:)<thresh);
    if ~isempty(active_times)
        if active_times(1) < (cutoff_time*12) % ignore cells who only turned on after 6+ hrs.
            metrics.off_times_nfkb(i) = active_times(end);
        end
    end    
end
metrics.off_times_nfkb = (metrics.off_times_nfkb-1)/12;
metrics.off_times_nfkb(metrics.off_times_nfkb<0) = 0;


%% NFkB OSCILLATION METRICS

% Calculate fourier distribution (via FFT) & power
Fs = 1/300;
depth = max(metrics.off_times_nfkb)*12;
NFFT = 2^nextpow2(depth); % Next power of 2 from chosen depth
aux.fft = zeros(size(metrics.time_series_nfkb,1),NFFT/2+1);
aux.freq = Fs/2*linspace(0,1,NFFT/2+1);
aux.power = zeros(size(aux.fft));


for i = 1:size(metrics.time_series_nfkb,1)
    if(metrics.off_times_nfkb(i)>0)
        y = metrics.time_series_nfkb(i,1:depth);
        off_frame = min([length(y), metrics.off_times_nfkb(i)*12+1+off_pad]); % (Pad w/ 1 extra hr of content)
        y(off_frame:end) = nan;
        y(isnan(y)) = [];
        y = y-nanmean(y);
        if ~isempty(y)
            Y = fft(y,NFFT)/length(y);
            aux.fft(i,:) = abs(Y(1:NFFT/2+1));
            aux.power(i,:) = abs(Y(1:NFFT/2+1).^2);
        end
    end
end

% Find the point of peak (secondary) power
metrics.peakfreq_nfkb = nan(size(aux.power,1),1);
for i =1:size(metrics.time_series_nfkb,1)
    [pks,locs] = globalpeaks(aux.power(i,:),2);
%     % (Code to check this "second-harmonic" thing)
%     if i<49
%     figure('Position',positionfig(220,100,[6,3])),
%     ha = tight_subplot(1,2);
%     plot(ha(1),1:100,metrics.time_series(i,1:100))
%     set(ha(1),'Ylim',[0 9],'XLim',[0 100],'Box','on')
%     hold(ha(2),'on')
%     plot(ha(2),freq, aux.power(i,:))
%     plot(ha(2), freq(locs),pks,'o')
%     hold(ha(2),'off')
%     set(ha(2),'XLim',[0 2],'Box','on')
%     end
    % Ensure we're not getting a totally spurious peak
    if min(pks) < (0.1*max(pks))
        locs(pks==min(pks)) = [];
    end
    if length(locs)>1
        idx = max(locs(1:2));
        metrics.peakfreq_nfkb(i) = 3600*aux.freq(idx);
    elseif ~isempty(locs)
         metrics.peakfreq_nfkb(i) = 3600*aux.freq(max([locs,3]));
    else
        metrics.peakfreq_nfkb(i) = 3600*aux.freq(1);
    end
end
%%
% Find total oscillatory content of particular cells (using thresholds from 0.35 to 0.7 hrs^(-1))
freq_thresh = aux.freq( (aux.freq >= (0.35/3600)) & (aux.freq <= (0.7/3600)));
metrics.oscfrac_nfkb = nan(size(aux.power,1),length(freq_thresh));
for j = 1:length(freq_thresh)
    for i =1:size(metrics.time_series_nfkb,1)
        metrics.oscfrac_nfkb(i,j) = nansum(aux.power(i,aux.freq >= freq_thresh(j))) /nansum(aux.power(i,:));
        if isnan(metrics.oscfrac_nfkb(i,j))
            metrics.oscfrac_nfkb(i,j) = 0;
        end
    end
end

%% NFkB METRICS OF AMPLITUDE AND TIMING
% 1st + 2nd peak time/amplitude
metrics.pk1_time_nfkb = nan(size(metrics.time_series_nfkb,1),1);
metrics.pk1_amp_nfkb =  nan(size(metrics.time_series_nfkb,1),1);
metrics.pk2_time_nfkb = nan(size(metrics.time_series_nfkb,1),1);
metrics.pk2_amp_nfkb =  nan(size(metrics.time_series_nfkb,1),1);
for i = 1:size(metrics.pk1_time_nfkb,1)    
    [pks, locs] = globalpeaks(metrics.time_series_nfkb(i,1:min([90,p.Results.MinLifetime])),5);
    % Supress any peaks that are within 6 frames of each other.
    [locs, order] = sort(locs,'ascend');
    pks = pks(order);
    while min(diff(locs))<6
        tmp = find(diff(locs)==min(diff(locs)),1,'first');
        tmp = tmp + (pks(tmp)>=pks(tmp+1));
        pks(tmp) = [];
        locs(tmp) = [];  
    end
    pks(locs<4) = [];
    locs(locs<4) = [];
    if ~isempty(locs)
        metrics.pk1_time_nfkb(i) = locs(1);
        metrics.pk1_amp_nfkb(i) = pks(1);
    end
    if length(locs)>1
        metrics.pk2_time_nfkb(i) = locs(2);
        metrics.pk2_amp_nfkb(i) = pks(2);
    end
end
metrics.pk1_time_nfkb = (metrics.pk1_time_nfkb-1)/12;
metrics.pk2_time_nfkb = (metrics.pk2_time_nfkb-1)/12;


%% NFkB METRICS OF DURATION
% Envelope width: maximum consecutive time above a threshold (envelope must begin within 1st 6 hrs)
smoothed2 = smoothrows(metrics.time_series_nfkb,5);
aux.thresholds = linspace(0, OnThreshNFkB*3, 40);
metrics.envelope_nfkb = zeros(size(metrics.time_series_nfkb,1),length(aux.thresholds));
for j = 1:length(aux.thresholds)
    thresholded = smoothed2>aux.thresholds(j);
    for i = 1:size(thresholded,1)
        curr = 1;
        idx_start = 1;
        while (curr<size(thresholded,2)) && (idx_start< (6*12))
            idx_start = find(thresholded(i,curr:end)==1,1,'first')+curr-1;
            if ~isempty(idx_start)
                idx_stop = find(thresholded(i,idx_start:end)==0,1,'first')+idx_start-1;
                if isempty(idx_stop)
                    idx_stop = find(~isnan(thresholded(i,:)),1,'last');
                end
                if (idx_stop-idx_start) > metrics.envelope_nfkb(i,j)
                    metrics.envelope_nfkb(i,j) = (idx_stop-idx_start);
                end
                curr = idx_stop;
            else
                break
            end
        end
    end
end
metrics.envelope_nfkb = metrics.envelope_nfkb/12;


% Number of frames above a given threshold
metrics.duration_nfkb = zeros(size(metrics.time_series_nfkb,1),length(aux.thresholds));
for i = 1:length(aux.thresholds)
    metrics.duration_nfkb(:,i) = nansum(smoothed>aux.thresholds(i),2)/12;
end
%% KTR METRICS
%% BASIC KTR METRICS: TIME SERIES, DERIVATIVE, INTEGRAL
% 1) basic time series. Interpolate over "normal" interval (12 frames per hr) if required
t = min(graph.t):1/12:max(graph.t);
if length(t)~=length(graph.t)
    metrics.time_series_ktr = nan(size(graph.var_ktr,1),length(t));
    for i = 1:size(graph.var_ktr,1)
        metrics.time_series_ktr(i,:) = interp1(graph.t,graph.var_ktr(i,:),t);
    end
else
    metrics.time_series_ktr = graph.var_ktr;
end

% 2) integrated activity
metrics.integrals_ktr = nan(size(metrics.time_series_ktr));
for i = 1:size(metrics.integrals_ktr,1)
    metrics.integrals_ktr(i,:) = cumtrapz(t,metrics.time_series_ktr(i,:));
end

% 3) differentiated activity - use central finite difference
smoothed = smoothrows(metrics.time_series_ktr,3);
metrics.derivatives_ktr = (smoothed(:,3:end) - smoothed(:,1:end-2))/(1/6);


%% TRIM EVERYBODY to a common length (of "good" sets, current minimum is roughly 21 hrs)
try
    metrics.time_series_ktr = metrics.time_series_ktr(:,1:p.Results.TrimFrame);
    metrics.integrals_ktr = metrics.integrals_ktr(:,1:p.Results.TrimFrame);
    metrics.derivatives_ktr = metrics.derivatives_ktr(:,1:(p.Results.TrimFrame-2));
    smoothed = smoothed(:,1:p.Results.TrimFrame);
    t = t(1:p.Results.TrimFrame);

catch me
    disp(['Note: vectors too short to cap @ ',num2str(p.Results.TrimFrame),' frames'])
end
%% MISC KTR METRICS

% 4) Integrals within one-hour windows (0-1, 1-2, 2-3) and three hour windows (0-3, 1-4, etc) of activity
max_hr = floor(max(t));
metrics.intwin1_ktr = nan(size(metrics.time_series_ktr,1),max_hr);
metrics.intwin3_ktr = nan(size(metrics.time_series_ktr,1),max_hr-2);
for i = 1:(max_hr)
    win = t>=(i-1) & t<(i);
    metrics.intwin1_ktr(:,i) = trapz(t(win),metrics.time_series_ktr(:,win),2);
    if i<= (max_hr-2)
        win = t>=(i-1) & t<(i+2);
        metrics.intwin3_ktr(:,i) = trapz(t(win),metrics.time_series_ktr(:,win),2);
    end
end

% MAX/MIN metrics
metrics.max_amplitude_ktr = nanmax(metrics.time_series_ktr,[],2);
metrics.max_integral_ktr = nanmax(metrics.integrals_ktr,[],2);
metrics.max_derivative_ktr = nanmax(metrics.derivatives_ktr,[],2);
metrics.min_derivative_ktr = nanmin(metrics.derivatives_ktr,[],2);

% ACTIVITY metrics: compute an off-time for all cells
metrics.off_times_ktr = zeros(size(smoothed,1),1);
inactive = [repmat(nanmin(smoothed(:,1:7),[],2),1,window_sz*2+1),smoothed(:,:),...
    repmat(nanmedian(smoothed(:,(end-window_sz:end)),2),1,window_sz*2)];
inactive = smoothrows(inactive<(OnThreshKTR),(window_sz*2));
frontcrop = round(window_sz*2*(1-thresh))+window_sz+1;
inactive = inactive(:,frontcrop:end);
inactive = inactive(:,1:size(smoothed,2));
inactive(isnan(smoothed)) = nan;

% Find the final time each cell was active
for i = 1:length(metrics.off_times_ktr)
    active_times = find(inactive(i,:)<thresh);
    if ~isempty(active_times)
        if active_times(1) < (cutoff_time*12) % ignore cells who only turned on after 6+ hrs.
            metrics.off_times_ktr(i) = active_times(end);
        end
    end    
end
metrics.off_times_ktr = (metrics.off_times_ktr-1)/12;
metrics.off_times_ktr(metrics.off_times_ktr<0) = 0;

%% KTR OSCILLATION METRICS
% Calculate fourier distribution (via FFT) & power
Fs = 1/300;
depth = max(metrics.off_times_ktr)*12;
NFFT = 2^nextpow2(depth); % Next power of 2 from chosen depth
aux.fft = zeros(size(metrics.time_series_ktr,1),NFFT/2+1);
aux.freq = Fs/2*linspace(0,1,NFFT/2+1);
aux.power = zeros(size(aux.fft));


for i = 1:size(metrics.time_series_ktr,1)
    if(metrics.off_times_ktr(i)>0)
        y = metrics.time_series_ktr(i,1:depth);
        off_frame = min([length(y), metrics.off_times_ktr(i)*12+1+off_pad]); % (Pad w/ 1 extra hr of content)
        y(off_frame:end) = nan;
        y(isnan(y)) = [];
        y = y-nanmean(y);
        if ~isempty(y)
            Y = fft(y,NFFT)/length(y);
            aux.fft(i,:) = abs(Y(1:NFFT/2+1));
            aux.power(i,:) = abs(Y(1:NFFT/2+1).^2);
        end
    end
end

% Find the point of peak (secondary) power
metrics.peakfreq_ktr = nan(size(aux.power,1),1);
for i =1:size(metrics.time_series_ktr,1)
    [pks,locs] = globalpeaks(aux.power(i,:),2);
%     % (Code to check this "second-harmonic" thing)
%     if i<49
%     figure('Position',positionfig(220,100,[6,3])),
%     ha = tight_subplot(1,2);
%     plot(ha(1),1:100,metrics.time_series(i,1:100))
%     set(ha(1),'Ylim',[0 9],'XLim',[0 100],'Box','on')
%     hold(ha(2),'on')
%     plot(ha(2),freq, aux.power(i,:))
%     plot(ha(2), freq(locs),pks,'o')
%     hold(ha(2),'off')
%     set(ha(2),'XLim',[0 2],'Box','on')
%     end
    % Ensure we're not getting a totally spurious peak
    if min(pks) < (0.1*max(pks))
        locs(pks==min(pks)) = [];
    end
    if length(locs)>1
        idx = max(locs(1:2));
        metrics.peakfreq_ktr(i) = 3600*aux.freq(idx);
    elseif ~isempty(locs)
         metrics.peakfreq_ktr(i) = 3600*aux.freq(max([locs,3]));
    else
        metrics.peakfreq_ktr(i) = 3600*aux.freq(1);
    end
end
%%
% Find total oscillatory content of particular cells (using thresholds from 0.35 to 0.7 hrs^(-1))
freq_thresh = aux.freq( (aux.freq >= (0.35/3600)) & (aux.freq <= (0.7/3600)));
metrics.oscfrac_ktr = nan(size(aux.power,1),length(freq_thresh));
for j = 1:length(freq_thresh)
    for i =1:size(metrics.time_series_ktr,1)
        metrics.oscfrac_ktr(i,j) = nansum(aux.power(i,aux.freq >= freq_thresh(j))) /nansum(aux.power(i,:));
        if isnan(metrics.oscfrac_ktr(i,j))
            metrics.oscfrac_ktr(i,j) = 0;
        end
    end
end

%% KTR METRICS OF AMPLITUDE AND TIMING
% 1st + 2nd peak time/amplitude
metrics.pk1_time_ktr = nan(size(metrics.time_series_ktr,1),1);
metrics.pk1_amp_ktr =  nan(size(metrics.time_series_ktr,1),1);
metrics.pk2_time_ktr = nan(size(metrics.time_series_ktr,1),1);
metrics.pk2_amp_ktr =  nan(size(metrics.time_series_ktr,1),1);
for i = 1:size(metrics.pk1_time_ktr,1)    
    [pks, locs] = globalpeaks(metrics.time_series_ktr(i,1:min([90,p.Results.MinLifetime])),5);
    % Supress any peaks that are within 6 frames of each other.
    [locs, order] = sort(locs,'ascend');
    pks = pks(order);
    while min(diff(locs))<6
        tmp = find(diff(locs)==min(diff(locs)),1,'first');
        tmp = tmp + (pks(tmp)>=pks(tmp+1));
        pks(tmp) = [];
        locs(tmp) = [];  
    end
    pks(locs<4) = [];
    locs(locs<4) = [];
    if ~isempty(locs)
        metrics.pk1_time_ktr(i) = locs(1);
        metrics.pk1_amp_ktr(i) = pks(1);
    end
    if length(locs)>1
        metrics.pk2_time_ktr(i) = locs(2);
        metrics.pk2_amp_ktr(i) = pks(2);
    end
end
metrics.pk1_time_ktr = (metrics.pk1_time_ktr-1)/12;
metrics.pk2_time_ktr = (metrics.pk2_time_ktr-1)/12;


%% KTR METRICS OF DURATION
% Envelope width: maximum consecutive time above a threshold (envelope must begin within 1st 6 hrs)
smoothed2 = smoothrows(metrics.time_series_ktr,5);
aux.thresholds = linspace(0, OnThreshKTR*3, 40);
metrics.envelope_ktr = zeros(size(metrics.time_series_ktr,1),length(aux.thresholds));
for j = 1:length(aux.thresholds)
    thresholded = smoothed2>aux.thresholds(j);
    for i = 1:size(thresholded,1)
        curr = 1;
        idx_start = 1;
        while (curr<size(thresholded,2)) && (idx_start< (6*12))
            idx_start = find(thresholded(i,curr:end)==1,1,'first')+curr-1;
            if ~isempty(idx_start)
                idx_stop = find(thresholded(i,idx_start:end)==0,1,'first')+idx_start-1;
                if isempty(idx_stop)
                    idx_stop = find(~isnan(thresholded(i,:)),1,'last');
                end
                if (idx_stop-idx_start) > metrics.envelope_ktr(i,j)
                    metrics.envelope_ktr(i,j) = (idx_stop-idx_start);
                end
                curr = idx_stop;
            else
                break
            end
        end
    end
end
metrics.envelope_ktr = metrics.envelope_ktr/12;


% Number of frames above a given threshold
metrics.duration_ktr = zeros(size(metrics.time_series_ktr,1),length(aux.thresholds));
for i = 1:length(aux.thresholds)
    metrics.duration_ktr(:,i) = nansum(smoothed>aux.thresholds(i),2)/12;
end
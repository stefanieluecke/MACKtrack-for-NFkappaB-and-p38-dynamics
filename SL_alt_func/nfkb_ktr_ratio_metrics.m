function [metrics,aux, graph, info, measure] = nfkb_ktr_ratio_metrics(id,varargin)
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
addParameter (p, 'OnThreshNFkB', 3, @isnumeric); %sigma threshold for determining responders
addParameter (p, 'OnThreshKTR', 3, @isnumeric); %sigma threshold for determining responders
addParameter(p, 'StartThreshNFkB', 14, @isnumeric); %max allowable starting threshhold (before baseline deduction)to filter out cells with pre-activated NFkB
addParameter(p,'StartThreshKTR',0.9, @isnumeric);
addParameter(p, 'MinSize', 90, @isnumeric); 
addParameter(p,'MinLifetime',109, @isnumeric);
addParameter(p,'TrimFrame',157, @isnumeric);
addParameter (p, 'GraphLimitsNFkB',[-0.25 7],@isnumeric);
addParameter (p, 'GraphLimitsKTR',[-0.02,0.35],@isnumeric);
expectedFlags = {'on','off'};
addParameter(p,'Verbose','off', @(x) any(validatestring(x,expectedFlags)))
valid_conv = @(x) assert(isnumeric(x)&&(x>=0)&&(length(x)==1),...
    'Convection correction parameter must be single integer >= 0');
addParameter(p,'ConvectionShift',1, valid_conv);
addParameter(p, 'StimulationTimePoint', 13, @isnumeric)
addParameter(p, 'FramesPerHour', 12, @isnumeric)
addParameter(p,'NFkBBaselineDeduction', 'on', @(x) any(validatestring(x,expectedFlags))) %option to turn off NFkB baseline deduction
addParameter(p, 'NFkBBackgroundAdjustment', 'on',@(x) any(validatestring(x,expectedFlags))) %option to turn off NFkB fluorescence distribution adjustment
addParameter(p,'NFkBBaselineAdjustment', 'on', @(x) any(validatestring(x,expectedFlags))) %option to turn off adjusment of NFkB trajectories with correction factor for fluorescence drop derived from Mock experiments


parse(p,id, varargin{:})

%% INITIALIZATION. Load and process data. Interpolate time series, calculate deriv/integral approximations
%{
%todo check if min life time is currently working fine
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
OnThreshNFkB = p.Results.OnThreshNFkB; % Minimum activity required for cell to register as 'on'
OnThreshKTR = p.Results.OnThreshKTR; % Minimum activity required for cell to register as 'on'
StimulationTimePoint = p.Results.StimulationTimePoint;
StartThreshNFkB = p.Results.StartThreshNFkB; 
StartThreshKTR = p.Results.StartThreshKTR; 
MinSize = p.Results.MinSize; 
MinLifetime = p.Results.MinLifetime; 
ConvectionShift = p.Results.ConvectionShift; 
FramesPerHour = p.Results.FramesPerHour;


%? double check all parameter passing through makes sense, esp baseline!
[graph, info, measure] = filter_nfkb_ktr_ratio(id,'MinLifetime',MinLifetime,...
                            'ConvectionShift',ConvectionShift, 'OnThreshNFkB',OnThreshNFkB,...
                            'OnThreshKTR',OnThreshKTR,'MinSize', MinSize,'StartThreshNFkB',...
                            StartThreshNFkB,'StartThreshKTR', StartThreshKTR, 'Verbose', p.Results.Verbose,...
                            'GraphLimitsNFkB', p.Results.GraphLimitsNFkB, 'GraphLimitsKTR', p.Results.GraphLimitsKTR, 'StimulationTimePoint', p.Results.StimulationTimePoint,...
                            'FramesPerHour', p.Results.FramesPerHour, 'NFkBBaselineDeduction', p.Results.NFkBBaselineDeduction, 'NFkBBackgroundAdjustment',p.Results.NFkBBackgroundAdjustment,'NFkBBaselineAdjustment', p.Results.NFkBBaselineAdjustment);
   
graph.var_nfkb = graph.var_nfkb(:,1:min(p.Results.TrimFrame, size(graph.var_nfkb,2))); %why is TrimFrame applied here and below?
graph.var_nfkb_no_base_ded = graph.var_nfkb_no_base_ded(:,1:min(p.Results.TrimFrame, size(graph.var_nfkb_no_base_ded,2)));
graph.var_ktr = graph.var_ktr(:,1:min(p.Results.TrimFrame, size(graph.var_ktr,2)));
graph.var_ktr_no_base_ded = graph.var_ktr_no_base_ded(:,1:min(p.Results.TrimFrame, size(graph.var_ktr_no_base_ded,2)));
graph.t = graph.t(1:size(graph.var_nfkb,2));

graph.opt_nfkb = maketicks(graph.t,info.GraphLimitsNFkB,0); %calls function to add tick labels to time frame (?only partially understand what that does)
graph.opt_ktr = maketicks(graph.t,info.GraphLimitsKTR,0); %calls function to add tick labels to time frame (?only partially understand what that does)
graph.opt_nfkb.Name = 'NF\kappaB activation'; 
graph.opt_ktr.Name = 'Kinase activation'; 

%do I need this?
if ~ismember ('OnThreshNFkB',p.UsingDefaults)
    OnThreshNFkB = info.OnThreshNFkB;
end
if ~ismember ('OnThreshKTR',p.UsingDefaults)
    OnThreshKTR = info.OnThreshKTR;
end

%% NFkB METRICS
%% BASIC NFkB METRICS: TIME SERIES, DERIVATIVE, INTEGRAL
% 1) basic time series. Interpolate over "normal" interval (12 frames per hr) if required
%    use baseline deducted NFkB, also generate non-baseline deducted NFkB here
t = min(graph.t):1/FramesPerHour:max(graph.t);
if length(t)~=length(graph.t)
    metrics.time_series_nfkb = nan(size(graph.var_nfkb,1),length(t));
    for i = 1:size(graph.var_nfkb,1)
        metrics.time_series_nfkb(i,:) = interp1(graph.t,graph.var_nfkb(i,:),t);
    end
    
    metrics.time_series_nfkb_no_base_ded = nan(size(graph.var_nfkb_no_base_ded,1),length(t));
    for i = 1:size(graph.var_nfkb_no_base_ded,1)
        metrics.time_series_nfkb_no_base_ded(i,:) = interp1(graph.t,graph.var_nfkb_no_base_ded(i,:),t);
    end
 else
    metrics.time_series_nfkb = graph.var_nfkb;
    metrics.time_series_nfkb_no_base_ded = graph.var_nfkb_no_base_ded;
end

%1.2) Include baseline calculated in filter function into metrics
metrics.baseline_nfkb = info.nfkb_baseline;

% 2) integrated activity
%includes below 0 activity 
% use baseline deducted time series
% start with stimulation time point

metrics.integrals_nfkb = cumtrapz(t(StimulationTimePoint:end),metrics.time_series_nfkb(:,StimulationTimePoint:end),2);

% 3) differentiated activity - use central finite difference
smoothed = smoothrows(metrics.time_series_nfkb,3);
metrics.derivatives_nfkb = (smoothed(:,3:end) - smoothed(:,1:end-2))/(1/6); %?divided by 1/6 because 2 tp, ie 10 min is 1/6 of an hour?


%% TRIM EVERYBODY to a common length (of "good" sets, current minimum is roughly 21 hrs)
try
    metrics.time_series_nfkb = metrics.time_series_nfkb(:,1:p.Results.TrimFrame);
    metrics.time_series_nfkb_no_base_ded = metrics.time_series_nfkb_no_base_ded(:,1:p.Results.TrimFrame);
    metrics.integrals_nfkb = metrics.integrals_nfkb(:,1:(p.Results.TrimFrame-StimulationTimePoint));
    metrics.derivatives_nfkb = metrics.derivatives_nfkb(:,1:(p.Results.TrimFrame-2));
    smoothed = smoothed(:,1:p.Results.TrimFrame);
    t = t(1:p.Results.TrimFrame);

catch me
    disp(['Note: vectors too short to cap @ ',num2str(p.Results.TrimFrame),' frames'])
end
%% MISC NFkB METRICS

% 4) Integrals within one-hour windows (0-1, 1-2, 2-3) after stimulation and three hour windows (0-3, 1-4, etc) of activity after stimulation
% use baseline deducted nfkb values
% automatically only includes time after stimulation (as long as t is adjusted to have 0 at StimulationTimePoint
max_hr = floor(max(t));
metrics.intwin1_nfkb = nan(size(metrics.time_series_nfkb,1),max_hr);
metrics.intwin3_nfkb = nan(size(metrics.time_series_nfkb,1),max_hr-2);
for i = 1:(max_hr)
    win = t>=(i-1) & t<(i); %eg first window goes from 0 to 1 h post stim 
    metrics.intwin1_nfkb(:,i) = trapz(t(win),metrics.time_series_nfkb(:,win),2);
    if i<= (max_hr-2)
        win = t>=(i-1) & t<(i+2); %eg first window goes from 0 to 3 h post stim 
        metrics.intwin3_nfkb(:,i) = trapz(t(win),metrics.time_series_nfkb(:,win),2);
    end
end

%% MAX/MIN metrics
%adjusted to include only time after stimulation
metrics.max_amplitude_nfkb = nanmax(metrics.time_series_nfkb(:,StimulationTimePoint:end),[],2);

%20200831 Testing shorter timeframe for max amplitude
%todo decide on timeframe
metrics.max_amplitude_4h_nfkb = nanmax(metrics.time_series_nfkb(:,StimulationTimePoint:48+StimulationTimePoint),[],2);

metrics.max_integral_nfkb = nanmax(metrics.integrals_nfkb,[],2);
metrics.max_derivative_nfkb = nanmax(metrics.derivatives_nfkb(:,StimulationTimePoint:end),[],2);
metrics.min_derivative_nfkb = nanmin(metrics.derivatives_nfkb(:,StimulationTimePoint:end),[],2);

%% ACTIVITY metrics:
%compute responder index using sigma threshold and off times
Wliml = StimulationTimePoint+1; %first/lower time point of window to check for activity
Wlimu = StimulationTimePoint + 48; %last/upper time point of window to check for activity, ie check in the first 4 hours after stimulation
blockLengthThresh = 5; %number of consecutive frames cell needs to pass activity threshold to be considered a responder
baseline_stdv_nfkb = nanstd(metrics.time_series_nfkb(:,1:StimulationTimePoint),0,2);

metrics.baseline_stdv_nfkb =baseline_stdv_nfkb;

NFkBBySigma = (smoothed(:,Wliml:Wlimu))./baseline_stdv_nfkb;
block_length_readout = zeros(size(NFkBBySigma,1),100);
block_length_readout(:,:) = 111;
for jj = 1:size(NFkBBySigma,1)
            thresh_idx = diff(NFkBBySigma(jj,:)>OnThreshNFkB,1,2);
            thresh_start = find(thresh_idx == 1);
            thresh_stop = find(thresh_idx == -1);
            if any(NFkBBySigma(jj,:)>OnThreshNFkB,2)
                if isempty(thresh_start) && isempty(thresh_stop)
                    block_length = Wlimu-Wliml;
                elseif isempty(thresh_start)
                    block_length = thresh_stop;
                elseif isempty(thresh_stop)
                    block_length = numel(thresh_idx)+1 - thresh_start;
                elseif ~isempty(thresh_start) && ~isempty(thresh_stop)
                    if (thresh_start(1)<thresh_stop(1)) && (thresh_stop(end)>thresh_start(end))
                           block_length = thresh_stop - thresh_start;
                    elseif thresh_start(1)<thresh_stop(1) && thresh_start(end)>thresh_stop(end)
                           block_length = [thresh_stop - thresh_start(1:end-1),numel(thresh_idx)+1 - thresh_start(end)];
                    elseif thresh_stop(1)<thresh_start(1) && thresh_stop(end)>thresh_start(end)
                           block_length = [thresh_stop(1),thresh_stop(2:end)-thresh_start];
                    elseif thresh_stop(1)<thresh_start(1) && thresh_start(end)>thresh_stop(end)
                           block_length = [thresh_stop(1),thresh_stop(2:end)-thresh_start(1:end-1), numel(thresh_idx)+1 - thresh_start(end)];
                    end                    
                end        
            else
                    block_length = 0;
            end 
            metrics.responder_index_nfkb(jj,1) = any(block_length>=blockLengthThresh, 2);
            block_length_readout(jj, 1:numel(block_length)) = block_length;
end

metrics.responders_fraction_nfkb = nnz(metrics.responder_index_nfkb)/numel(metrics.responder_index_nfkb);
%%
% Determination of off_times
% using smoothed values and stdv of non-smoothed basline 
smoothed_by_sigma = smoothed./baseline_stdv_nfkb;
on_array = zeros(size(smoothed(:,Wliml:end)));
metrics.off_times_nfkb = zeros(size(smoothed,1),1);
for ii = 1:size(smoothed_by_sigma,1)
    n = 0;
    for jj = Wliml:size(smoothed_by_sigma,2)
        if smoothed_by_sigma(ii,jj)> OnThreshNFkB
            n = n+1;
        else
        n = 0;
        end
        on_array(ii,(jj-Wliml+1)) = n;
    end
    if find(on_array(ii,:)==1, 1)> Wlimu %ignore cells activating for first time after expected activity window
       metrics.off_times_nfkb(ii) = 0;
    else
        if ~isempty(find(on_array(ii,:)>= blockLengthThresh, 1, 'last'))
           metrics.off_times_nfkb(ii) = find(on_array(ii,:)>=blockLengthThresh, 1, 'last');
        else
            metrics.off_times_nfkb(ii) = 0;
        end
    end
end
%metrics.off_times_nfkb = (metrics.off_times_nfkb-StimulationTimePoint)/FramesPerHour;
metrics.off_times_nfkb = metrics.off_times_nfkb/FramesPerHour; %adjustment for stimulation time point not necessary here, because we only start counting frames after stimulation (winl)
metrics.off_times_nfkb(metrics.off_times_nfkb<0) = 0;

%% NFkB OSCILLATION METRICS

% Calculate fourier distribution (via FFT) & power
%todo check if off_pad is needed
off_pad = 12; % Signal time added to trajectory in  FFT calculation (keeps transients from being recorded as osc.)
Fs = 1/300;
depth = max(metrics.off_times_nfkb)*FramesPerHour;
NFFT = 2^nextpow2(depth); % Next power of 2 from chosen depth, can be used to pad input to fft
aux.nfkb.fft = zeros(size(metrics.time_series_nfkb,1),NFFT/2+1);
aux.nfkb.freq = Fs/2*linspace(0,1,NFFT/2+1);
aux.nfkb.power = zeros(size(aux.nfkb.fft));
aux.nfkb.power_norm = zeros(size(aux.nfkb.fft));


for i = 1:size(metrics.time_series_nfkb,1)
    if(metrics.off_times_nfkb(i)>0)
        %adjusted to start with stimulation
        y = metrics.time_series_nfkb(i,StimulationTimePoint:(depth+StimulationTimePoint));
        off_frame = min([length(y), metrics.off_times_nfkb(i)*FramesPerHour+1+off_pad]); % (Pad w/ 1 extra hr of content)
        y(off_frame:end) = nan;
        y(isnan(y)) = [];
        y = y-nanmean(y);
        if ~isempty(y)
            %SL 20200824: remove division by length(y), to be able to use Plancherel's theorem correctly below
           Y = fft(y,NFFT);
%            Y = fft(y,NFFT)/length(y); %20200902 SL What was the division by y needed/used for?
            aux.nfkb.fft(i,:) = abs(Y(1:NFFT/2+1));
            aux.nfkb.power(i,:) = abs(Y(1:NFFT/2+1).^2);
             %SL20200820 Normalization (Plancherel's theorem) to make power peaks comaparable
            aux.nfkb.power_norm(i,:) = abs(Y(1:NFFT/2+1).^2)/sum(abs(y).^2);

        end
    end
end

% Original Peakfreq metric using aux.power

% Find the point of peak (secondary) power
metrics.peakfreq_nfkb = nan(size(aux.nfkb.power,1),1);
for i =1:size(metrics.time_series_nfkb,1)
    [pks,locs] = globalpeaks(aux.nfkb.power(i,:),2);
%     % (Code to check this "second-harmonic" thing)
%     if i<49
%     figure('Position',positionfig(220,100,[6,3])),
%     ha = tight_subplot(1,2);
%     plot(ha(1),1:100,metrics.time_series(i,1:100))
%     set(ha(1),'Ylim',[0 9],'XLim',[0 100],'Box','on')
%     hold(ha(2),'on')
%     plot(ha(2),freq, aux.nfkb.power(i,:))
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
        metrics.peakfreq_nfkb(i) = 3600*aux.nfkb.freq(idx);
    elseif ~isempty(locs)
         metrics.peakfreq_nfkb(i) = 3600*aux.nfkb.freq(max([locs,3]));
    else
        metrics.peakfreq_nfkb(i) = 3600*aux.nfkb.freq(1);
    end
end

% New Peakfreq metric using normalized aux.power_norm

% Find the point of peak (secondary) power
metrics.peakfreq_norm_nfkb = nan(size(aux.nfkb.power,1),1);
for i =1:size(metrics.time_series_nfkb,1)
    [pks,locs] = globalpeaks(aux.nfkb.power_norm(i,:),2);
%     % (Code to check this "second-harmonic" thing)
%     if i<49
%     figure('Position',positionfig(220,100,[6,3])),
%     ha = tight_subplot(1,2);
%     plot(ha(1),1:100,metrics.time_series(i,1:100))
%     set(ha(1),'Ylim',[0 9],'XLim',[0 100],'Box','on')
%     hold(ha(2),'on')
%     plot(ha(2),freq, aux.nfkb.power(i,:))
%     plot(ha(2), freq(locs),pks,'o')
%     hold(ha(2),'off')
%     set(ha(2),'XLim',[0 2],'Box','on')
%     end
    % Ensure we're not getting a totally spurious peak
    if min(pks) < (0.1*max(pks))
        locs(pks==min(pks)) = [];
    elseif min(pks)< 6.5
        locs(pks==min(pks))=[];
    end
    
    if max(pks)< 6.5
        locs=[];
    end
    
    if length(locs)>1
        idx = max(locs(1:2));
        metrics.peakfreq_norm_nfkb(i) = 3600*aux.nfkb.freq(idx);
    elseif ~isempty(locs)
         metrics.peakfreq_norm_nfkb(i) = 3600*aux.nfkb.freq(max([locs,3]));
    else
        metrics.peakfreq_norm_nfkb(i) = 3600*aux.nfkb.freq(1);
    end
end

%%
% Find total oscillatory content of particular cells (using thresholds from 0.35 to 0.7 hrs^(-1))
freq_thresh = aux.nfkb.freq((aux.nfkb.freq >= (0.35/3600)) & (aux.nfkb.freq <= (0.7/3600)));
%uisng non-normalized power
metrics.oscfrac_nfkb = nan(size(aux.nfkb.power,1),length(freq_thresh));
for j = 1:length(freq_thresh)
    for i =1:size(metrics.time_series_nfkb,1)
        metrics.oscfrac_nfkb(i,j) = nansum(aux.nfkb.power(i,aux.nfkb.freq >= freq_thresh(j))) /nansum(aux.nfkb.power(i,:));
        if isnan(metrics.oscfrac_nfkb(i,j))
            metrics.oscfrac_nfkb(i,j) = 0;
        end
    end
end

%uisng normalized power
metrics.oscfrac_norm_nfkb = nan(size(aux.nfkb.power_norm,1),length(freq_thresh));
for j = 1:length(freq_thresh)
    for i =1:size(metrics.time_series_nfkb,1)
        metrics.oscfrac_norm_nfkb(i,j) = nansum(aux.nfkb.power_norm(i,((aux.nfkb.freq >= freq_thresh(j))&(aux.nfkb.power_norm(i,:)>= 6.5)))) /nansum(aux.nfkb.power_norm(i,aux.nfkb.power_norm(i,:)>= 6.5));
%        metrics.oscfrac_norm_nfkb(i,j) = nansum(aux.nfkb.power_norm(i,aux.nfkb.freq >= freq_thresh(j))) /nansum(aux.nfkb.power_norm(i,:));
        if isnan(metrics.oscfrac_norm_nfkb(i,j))
            metrics.oscfrac_norm_nfkb(i,j) = 0;
        end
    end
end

%% NFkB METRICS OF AMPLITUDE AND TIMING

% 1st + 2nd peak time/amplitude/prominence/width/height
pk_feats = {'pk1_amp_nfkb', 'pk1_time_nfkb', 'pk1_width_nfkb', 'pk1_prom_nfkb', 'pk1_height_nfkb',...
        'pk2_amp_nfkb', 'pk2_time_nfkb', 'pk2_width_nfkb', 'pk2_prom_nfkb', 'pk2_height_nfkb'};
for i=1:length(pk_feats)
    metrics.(pk_feats{i}) = nan(size(metrics.time_series_nfkb,1),1);
end

for i = 1:size(metrics.pk1_time_nfkb,1)
    
    %todo testing if smoothing is helpful
    %smoothing of trajectories
 %   metrics.time_series_smoothed_nfkb = smoothrows(metrics.time_series_nfkb,3);
    
  %  [pks, locs] =  globalpeaks(metrics.time_series_nfkb(i,1:min([90,p.Results.MinLifetime])),5);
    %todo check if I want to stick with 90 TPs (7.5 h, ie 6.5 h + before
    %stimulation), test shorter times
  %  %20200617 test to include more peaks because of more filtering
%    [pks_nfkb, locs_nfkb, width_nfkb, prom_nfkb, heights_nfkb] = globalpeaks(metrics.time_series_smoothed_nfkb(i,1:min([48+StimulationTimePoint,p.Results.MinLifetime])),5);
     [pks_nfkb, locs_nfkb, width_nfkb, prom_nfkb, heights_nfkb] = globalpeaks(metrics.time_series_nfkb(i,1:min([48+StimulationTimePoint,p.Results.MinLifetime])),10);
%    [pks_nfkb, locs_nfkb, width_nfkb, prom_nfkb, heights_nfkb] = globalpeaks(metrics.time_series_nfkb(i,1:min([96+StimulationTimePoint,p.Results.MinLifetime])),5);
%    [pks_nfkb, locs_nfkb, width_nfkb, prom_nfkb, heights_nfkb] = globalpeaks(metrics.time_series_nfkb(i,1:min([90,p.Results.MinLifetime])),5);
    % Supress any peaks that are within 6 frames of each other.
    [locs_nfkb, order_nfkb] = sort(locs_nfkb,'ascend');
    pks_nfkb = pks_nfkb(order_nfkb);width_nfkb = width_nfkb(order_nfkb); prom_nfkb = prom_nfkb(order_nfkb); heights_nfkb = heights_nfkb(order_nfkb);
   
    while min(diff(locs_nfkb))<6
        tmp = find(diff(locs_nfkb)==min(diff(locs_nfkb)),1,'first');
        tmp = tmp + (pks_nfkb(tmp)>=pks_nfkb(tmp+1));
        pks_nfkb(tmp) = []; locs_nfkb(tmp) = []; width_nfkb(tmp) = []; prom_nfkb(tmp) = []; heights_nfkb(tmp) = [];
    end
    
    pks_nfkb(locs_nfkb<(StimulationTimePoint + 1)) = [];
    width_nfkb(locs_nfkb<(StimulationTimePoint + 1)) = [];
    prom_nfkb(locs_nfkb<(StimulationTimePoint + 1)) = [];
    heights_nfkb(locs_nfkb<(StimulationTimePoint + 1)) = [];
    locs_nfkb(locs_nfkb<(StimulationTimePoint + 1)) = [];
%    pks(locs<(StimulationTimePoint + 3)) = []; %this seems to make some early peaks not detected --> remove extra padding
%    locs(locs<(StimulationTimePoint + 3)) = [];

% 20200617 Testing: Filter to remove peaks with amplitudes below 0
%todo test this
    locs_nfkb(pks_nfkb<= 0) = [];
    width_nfkb(pks_nfkb<= 0) = [];
    prom_nfkb(pks_nfkb<= 0) = [];
    heights_nfkb(pks_nfkb<= 0) = [];
    pks_nfkb(pks_nfkb<= 0) = []; %pks needs to be filtered after others
%20200826 SL Testing: Filter to remove peaks with short peak prominence based on Stdv of baseline
 %   locs_nfkb(prom_nfkb< 2*baseline_stdv_nfkb(i)) = [];
 %   width_nfkb(prom_nfkb< 2*baseline_stdv_nfkb(i)) = [];
  %  pks_nfkb(prom_nfkb< 2*baseline_stdv_nfkb(i)) = []; 
   % heights_nfkb(prom_nfkb< 2*baseline_stdv_nfkb(i)) = [];
    %prom_nfkb(prom_nfkb< 2*baseline_stdv_nfkb(i)) = [];%prom pks needs to be filtered after others

    locs_nfkb(heights_nfkb< 2*baseline_stdv_nfkb(i)) = [];
    width_nfkb(heights_nfkb< 2*baseline_stdv_nfkb(i)) = [];
    pks_nfkb(heights_nfkb< 2*baseline_stdv_nfkb(i)) = []; 
    prom_nfkb(heights_nfkb< 2*baseline_stdv_nfkb(i)) = [];
    heights_nfkb(heights_nfkb< 2*baseline_stdv_nfkb(i)) = [];%heights pks needs to be filtered after others
 
%
    %20200923 SL Testing: Filter to remove peaks that are too narrow (based on 'halfheight' width determination in globablpeaks)
    locs_nfkb(width_nfkb< 2) = [];
    pks_nfkb(width_nfkb< 2) = []; 
    prom_nfkb(width_nfkb< 2) = [];
    heights_nfkb(width_nfkb< 2)= []; 
    width_nfkb(width_nfkb< 2)= [];%width of pks needs to be filtered last 
  %}  
    
   if ~isempty(locs_nfkb)
        metrics.pk1_time_nfkb(i) = locs_nfkb(1);
        metrics.pk1_amp_nfkb(i) = pks_nfkb(1);
        metrics.pk1_width_nfkb(i) = width_nfkb(1);
        metrics.pk1_prom_nfkb(i) = prom_nfkb(1);
        metrics.pk1_height_nfkb(i) = heights_nfkb(1);

    end
    if length(locs_nfkb)>1
        metrics.pk2_time_nfkb(i) = locs_nfkb(2);
        metrics.pk2_amp_nfkb(i) = pks_nfkb(2);
        metrics.pk2_width_nfkb(i) = width_nfkb(2);
        metrics.pk2_prom_nfkb(i) = prom_nfkb(2);
        metrics.pk2_height_nfkb(i) = heights_nfkb(2);

    end
end
metrics.pk1_time_nfkb = (metrics.pk1_time_nfkb-StimulationTimePoint)/FramesPerHour;
metrics.pk2_time_nfkb = (metrics.pk2_time_nfkb-StimulationTimePoint)/FramesPerHour;
%metrics.pk1_time_nfkb = (metrics.pk1_time_nfkb-1)/12;
%metrics.pk2_time_nfkb = (metrics.pk2_time_nfkb-1)/12;


%% NFkB METRICS OF DURATION
% Envelope width: maximum consecutive time above a threshold (envelope must begin within 1st 6 hrs)
smoothed2 = smoothrows(metrics.time_series_nfkb,5);
%aux.nfkb.thresholds = linspace(0, OnThreshNFkB*3, 40);
upperThresh = 7.1;
aux.nfkb.thresholds = linspace(0, upperThresh, 25);
metrics.envelope_nfkb = zeros(size(metrics.time_series_nfkb,1),length(aux.nfkb.thresholds));
for j = 1:length(aux.nfkb.thresholds)
    thresholded = smoothed2(:,StimulationTimePoint:end)>aux.nfkb.thresholds(j);
    for i = 1:size(thresholded,1)
        curr = 1;
        idx_start = 1;
        while (curr<size(thresholded,2)) && (idx_start< (6*FramesPerHour))
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
metrics.envelope_nfkb = metrics.envelope_nfkb/FramesPerHour;


% Number of frames above a given threshold
metrics.duration_nfkb = zeros(size(metrics.time_series_nfkb,1),length(aux.nfkb.thresholds));
for i = 1:length(aux.nfkb.thresholds)
    metrics.duration_nfkb(:,i) = nansum(smoothed(:,StimulationTimePoint:end)>aux.nfkb.thresholds(i),2)/FramesPerHour;
end

%% KTR METRICS
%% BASIC KTR METRICS: TIME SERIES, Baseline DERIVATIVE, INTEGRAL
%
% 1.1) basic time series. Interpolate over "normal" interval (12 frames per hr) if required
% generate time series both with non-baseline deducted and baseline-deducted
t = min(graph.t):1/FramesPerHour:max(graph.t);
if length(t)~=length(graph.t)
    metrics.time_series_ktr = nan(size(graph.var_ktr,1),length(t));
    for i = 1:size(graph.var_ktr,1)
        metrics.time_series_ktr(i,:) = interp1(graph.t,graph.var_ktr(i,:),t);
    end
    
       metrics.time_series_ktr_no_base_ded = nan(size(graph.var_ktr_no_base_ded,1),length(t));
    for i = 1:size(graph.var_ktr_no_base_ded,1)
        metrics.time_series_ktr_no_base_ded(i,:) = interp1(graph.t,graph.var_ktr_no_base_ded(i,:),t);
    end
    
else
    metrics.time_series_ktr = graph.var_ktr;
    metrics.time_series_ktr_no_base_ded = graph.var_ktr_no_base_ded;
end
%1.2) Include baseline calculated in filter function into metrics
metrics.baseline_ktr = info.ktr_baseline;

% 2) integrated activity
%use baseline deducted KTR ratios for this
%integration incl values below 0
    metrics.integrals_ktr = cumtrapz(t(StimulationTimePoint:end),metrics.time_series_ktr(:, StimulationTimePoint:end),2);

% 3) differentiated activity - use central finite difference
% use baseline deducted KTR ratios for this
smoothed = smoothrows(metrics.time_series_ktr,3);
metrics.derivatives_ktr = (smoothed(:,3:end) - smoothed(:,1:end-2))/(1/6);


%% TRIM EVERYBODY to a common length (of "good" sets, current minimum is roughly 21 hrs)
try
    metrics.time_series_ktr = metrics.time_series_ktr(:,1:p.Results.TrimFrame);
    metrics.time_series_ktr_no_base_ded = metrics.time_series_ktr_no_base_ded(:,1:p.Results.TrimFrame);
    metrics.integrals_ktr = metrics.integrals_ktr(:,1:(p.Results.TrimFrame-StimulationTimePoint));
    metrics.derivatives_ktr = metrics.derivatives_ktr(:,1:(p.Results.TrimFrame-2));
    smoothed = smoothed(:,1:p.Results.TrimFrame);
    t = t(1:p.Results.TrimFrame);

catch me
    disp(['Note: vectors too short to cap @ ',num2str(p.Results.TrimFrame),' frames'])
end
%% MISC KTR METRICS

% 4) Integrals within one-hour windows (0-1, 1-2, 2-3) after stimulation and three hour windows (0-3, 1-4, etc) of activity after stimulation
% use baseline deducted ktr values
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

%% MAX/MIN metrics
metrics.max_amplitude_ktr = nanmax(metrics.time_series_ktr(:,StimulationTimePoint:end),[],2);

%20200831 Testing shorter timeframe for max amplitude
%todo decide on timeframe
metrics.max_amplitude_4h_ktr = nanmax(metrics.time_series_ktr(:,StimulationTimePoint:48+StimulationTimePoint),[],2);

metrics.max_integral_ktr = nanmax(metrics.integrals_ktr,[],2);
metrics.max_derivative_ktr = nanmax(metrics.derivatives_ktr(:,StimulationTimePoint:end),[],2);
metrics.min_derivative_ktr = nanmin(metrics.derivatives_ktr(:,StimulationTimePoint:end),[],2);

%% ACTIVITY metrics: compute an off-time for all cells

%compute responder index using sigma threshold and off times
Wliml = StimulationTimePoint+1; %first/lower time point of window to check for activity
Wlimu = StimulationTimePoint + 36; %last/upper time point of window to check for activity, ie check in the first 3 hours after stimulation
blockLengthThresh = 3; %number of consecutive frames cell needs to pass activity threshold to be considered a responder
baseline_stdv_ktr = nanstd(metrics.time_series_ktr(:,1:StimulationTimePoint),0,2);

metrics.baseline_stdv_ktr =baseline_stdv_ktr;

ktrBySigma = (smoothed(:,Wliml:Wlimu))./baseline_stdv_ktr;
block_length_readout = zeros(size(ktrBySigma,1),100);
block_length_readout(:,:) = 111;
for jj = 1:size(ktrBySigma,1)
            thresh_idx = diff(ktrBySigma(jj,:)>OnThreshKTR,1,2);
            thresh_start = find(thresh_idx == 1);
            thresh_stop = find(thresh_idx == -1);
            if any(ktrBySigma(jj,:)>OnThreshKTR,2)
                if isempty(thresh_start) && isempty(thresh_stop)
                    block_length = Wlimu-Wliml;
                elseif isempty(thresh_start)
                    block_length = thresh_stop;
                elseif isempty(thresh_stop)
                    block_length = numel(thresh_idx)+1 - thresh_start;
                elseif ~isempty(thresh_start) && ~isempty(thresh_stop)
                    if (thresh_start(1)<thresh_stop(1)) && (thresh_stop(end)>thresh_start(end))
                           block_length = thresh_stop - thresh_start;
                    elseif thresh_start(1)<thresh_stop(1) && thresh_start(end)>thresh_stop(end)
                           block_length = [thresh_stop - thresh_start(1:end-1),numel(thresh_idx)+1 - thresh_start(end)];
                    elseif thresh_stop(1)<thresh_start(1) && thresh_stop(end)>thresh_start(end)
                           block_length = [thresh_stop(1),thresh_stop(2:end)-thresh_start];
                    elseif thresh_stop(1)<thresh_start(1) && thresh_start(end)>thresh_stop(end)
                           block_length = [thresh_stop(1),thresh_stop(2:end)-thresh_start(1:end-1), numel(thresh_idx)+1 - thresh_start(end)];
                    end                    
                end        
            else
                    block_length = 0;
            end 
            metrics.responder_index_ktr(jj,1) = any(block_length>=blockLengthThresh, 2);
            block_length_readout(jj, 1:numel(block_length)) = block_length;
end

metrics.responders_fraction_ktr = nnz(metrics.responder_index_ktr)/numel(metrics.responder_index_ktr);

%%
%Off times

%using smoothed data and stdv of non-smoothedbasline

smoothed_by_sigma = smoothed./baseline_stdv_ktr;
on_array = zeros(size(smoothed(:,Wliml:end)));
metrics.off_times_ktr = zeros(size(smoothed,1),1);
for ii = 1:size(smoothed_by_sigma,1)
    n = 0;
    for jj = Wliml:size(smoothed_by_sigma,2)
        if smoothed_by_sigma(ii,jj)> OnThreshKTR
            n = n+1;
        else
        n = 0;
        end
        on_array(ii,(jj-Wliml+1)) = n;
    end
    if find(on_array(ii,:)==1, 1)> Wlimu %ignore cells activating for first time after expected activity window
       metrics.off_times_ktr(ii) = 0;
    else
        if ~isempty(find(on_array(ii,:)>= blockLengthThresh, 1, 'last'))
           metrics.off_times_ktr(ii) = find(on_array(ii,:)>=blockLengthThresh, 1, 'last');
        else
            metrics.off_times_ktr(ii) = 0;
        end
    end
end
metrics.off_times_ktr = metrics.off_times_ktr/FramesPerHour;
metrics.off_times_ktr(metrics.off_times_ktr<0) = 0;

%% KTR OSCILLATION METRICS
% Calculate fourier distribution (via FFT) & power
%todo check of off_pad is needed with my new offtimes, etc
off_pad = FramesPerHour; % Signal time added to trajectory in  FFT calculation (keeps transients from being recorded as osc.)
Fs = 1/300;
depth = max(metrics.off_times_ktr)*FramesPerHour;
NFFT = 2^nextpow2(depth); % Next power of 2 from chosen depth
aux.ktr.fft = zeros(size(metrics.time_series_ktr,1),NFFT/2+1);
aux.ktr.freq = Fs/2*linspace(0,1,NFFT/2+1);
aux.ktr.power = zeros(size(aux.ktr.fft));
aux.ktr.power_norm = zeros(size(aux.ktr.fft));
for i = 1:size(metrics.time_series_ktr,1)
    if(metrics.off_times_ktr(i)>0) %only calculate peakfreq if offtime is not 0
        y = metrics.time_series_ktr(i,StimulationTimePoint:(depth+StimulationTimePoint));
%       y = metrics.time_series_ktr(i,1:depth);
        %todo see if I need to increase padding due to very strict offtime determination
        %todo why don't we use the entire trajectory?
        off_frame = min([length(y), metrics.off_times_ktr(i)*FramesPerHour+1+off_pad]); % (Pad w/ 1 extra hr of content)
        y(off_frame:end) = nan;
        y(isnan(y)) = [];
        y = y-nanmean(y);
        if ~isempty(y)
            %SL 20200824: remove division by length(y), to be able to use Plancherel's theorem correctly below
            Y = fft(y,NFFT);
%           Y = fft(y,NFFT)/length(y);
            aux.ktr.fft(i,:) = abs(Y(1:NFFT/2+1));
            aux.ktr.power(i,:) = abs(Y(1:NFFT/2+1).^2);
            %SL20200820 Normalization (Plancherel's theorem) to make power peaks comaparable
           aux.ktr.power_norm(i,:) = abs(Y(1:NFFT/2+1).^2)/sum((abs(y)).^2);
        end
    end
end

%todo decide on normalized or non-normalized method and remove the rest from the code
% Original Peakfreq metric using aux.power
% Find the point of peak (secondary) power
metrics.peakfreq_ktr = nan(size(aux.ktr.power,1),1);
for i =1:size(metrics.time_series_ktr,1)
    [pks,locs] = globalpeaks(aux.ktr.power(i,:),2);
%     % (Code to check this "second-harmonic" thing)
%     if i<49
%     figure('Position',positionfig(220,100,[6,3])),
%     ha = tight_subplot(1,2);
%     plot(ha(1),1:100,metrics.time_series(i,1:100))
%     set(ha(1),'Ylim',[0 9],'XLim',[0 100],'Box','on')
%     hold(ha(2),'on')
%     plot(ha(2),freq, aux.ktr.power(i,:))
%     plot(ha(2), freq(locs),pks,'o')
%     hold(ha(2),'off')
%     set(ha(2),'XLim',[0 2],'Box','on')
%     end
    % Ensure we're not getting a totally spurious peak
    
    %todo check what to do/how to deal with a number of similarly sized
    %peaks in fft that represent 'oscillating' noise!
    
    %if the two found peaks are not obviously different heights, both peaks
    %kept
    if min(pks) < (0.1*max(pks))
        locs(pks==min(pks)) = [];
    end
    %if both peaks kept, use higher frequency peak as peakfreq
    if length(locs)>1
        idx = max(locs(1:2));
        metrics.peakfreq_ktr(i) = 3600*aux.ktr.freq(idx);
    %if locs has 1 value, it takes either frequency corresponding to peak
    %or 3rd frequency (close to 0), whichever is larger
    elseif ~isempty(locs)
         metrics.peakfreq_ktr(i) = 3600*aux.ktr.freq(max([locs,3]));
    %if no peaks, peakfreq is 0
    else
        metrics.peakfreq_ktr(i) = 3600*aux.ktr.freq(1);
    end
end


% New Peakfreq metric using normalized aux.power_norm
% Find the point of peak (secondary) power
metrics.peakfreq_norm_ktr = nan(size(aux.ktr.power_norm,1),1);
for i =1:size(metrics.time_series_ktr,1)
    [pks,locs] = globalpeaks(aux.ktr.power_norm(i,:),2);
%     % (Code to check this "second-harmonic" thing)
%     if i<49
%     figure('Position',positionfig(220,100,[6,3])),
%     ha = tight_subplot(1,2);
%     plot(ha(1),1:100,metrics.time_series(i,1:100))
%     set(ha(1),'Ylim',[0 9],'XLim',[0 100],'Box','on')
%     hold(ha(2),'on')
%     plot(ha(2),freq, aux.ktr.power(i,:))
%     plot(ha(2), freq(locs),pks,'o')
%     hold(ha(2),'off')
%     set(ha(2),'XLim',[0 2],'Box','on')
%     end
    % Ensure we're not getting a totally spurious peak
    
    %todo check what to do/how to deal with a number of similarly sized
    %peaks in fft that represent 'oscillating' noise!
    
    %if the two found peaks are not obviously different heights, both peaks
    %kept
    if min(pks) < (0.1*max(pks))
        locs(pks==min(pks)) = [];
    elseif min(pks)< 6.5
        locs(pks==min(pks))=[];
    end
    
    if max(pks)< 6.5
        %locs(pks==max(pks))=[];
        locs=[];
    end

    %if both peaks kept, use higher frequency peak as peakfreq
    if length(locs)>1
        idx = max(locs(1:2));
        metrics.peakfreq_norm_ktr(i) = 3600*aux.ktr.freq(idx);
    %if locs has 1 value, it takes either frequency corresponding to peak
    %or 3rd frequency (close to 0), whichever is larger, to avoid mixing
    %the low peakfreq determined here becoming indistinguishable from the
    %non-responders with peakfreq = 0

    elseif ~isempty(locs)
         metrics.peakfreq_norm_ktr(i) = 3600*aux.ktr.freq(max([locs,3]));
    %if no peaks, peakfreq is 0
    else
        metrics.peakfreq_norm_ktr(i) = 3600*aux.ktr.freq(1);
    end
end


%%
% Find total oscillatory content of particular cells (using thresholds from 0.38 to 4 hrs^(-1))

%20200819 Adjust thresholds for expected frequency range in KTR
%20200824 After adding normalization and cutoff value to power, go back to previous freq range
%freq_thresh = aux.ktr.freq( (aux.ktr.freq >= (0.38/3600)) & (aux.ktr.freq <= (4/3600)));
freq_thresh = aux.ktr.freq( (aux.ktr.freq >= (0.35/3600)) & (aux.ktr.freq <= (0.7/3600)));

%using non-normalized power
metrics.oscfrac_ktr = nan(size(aux.ktr.power,1),length(freq_thresh));
for j = 1:length(freq_thresh)
    for i =1:size(metrics.time_series_ktr,1)
        metrics.oscfrac_ktr(i,j) = nansum(aux.ktr.power(i,aux.ktr.freq >= freq_thresh(j))) /nansum(aux.ktr.power(i,:));
        if isnan(metrics.oscfrac_ktr(i,j))
            metrics.oscfrac_ktr(i,j) = 0;
        end
    end
end


%using normalized power
metrics.oscfrac_norm_ktr = nan(size(aux.ktr.power_norm,1),length(freq_thresh));
for j = 1:length(freq_thresh)
    for i =1:size(metrics.time_series_ktr,1)
        metrics.oscfrac_norm_ktr(i,j) = nansum(aux.ktr.power_norm(i,((aux.ktr.freq >= freq_thresh(j))&(aux.ktr.power_norm(i,:)>= 6.5)))) /nansum(aux.ktr.power_norm(i,aux.ktr.power_norm(i,:)>= 6.5));
%        metrics.oscfrac_norm_ktr(i,j) = nansum(aux.ktr.power_norm(i,aux.ktr.freq >= freq_thresh(j))) /nansum(aux.ktr.power_norm(i,:));
        if isnan(metrics.oscfrac_norm_ktr(i,j))
            metrics.oscfrac_norm_ktr(i,j) = 0;
        end
    end
end
%}

%% KTR METRICS OF AMPLITUDE AND TIMING
% 1st + 2nd peak time/amplitude/prominence/width/height
pk_feats = {'pk1_amp_ktr', 'pk1_time_ktr', 'pk1_width_ktr', 'pk1_prom_ktr', 'pk1_height_ktr',...
        'pk2_amp_ktr', 'pk2_time_ktr', 'pk2_width_ktr', 'pk2_prom_ktr', 'pk2_height_ktr'};
for i=1:length(pk_feats)
    metrics.(pk_feats{i}) = nan(size(metrics.time_series_ktr,1),1);
end

for i = 1:size(metrics.pk1_time_ktr,1)    

    %todo testing if smoothing is helpful
    %smoothing of trajectories
 %   metrics.time_series_smoothed_ktr = smoothrows(metrics.time_series_ktr,3);

    
    %  globalpeaks(metrics.time_series_ktr(i,1:min([90,p.Results.MinLifetime])),5);
    %todo check if I want to stick with 90 TPs (7.5 h, ie 6.5 h + before stimulation), test shorter time frames
  %  %20200617 test to include more peaks because of more filtering
%    [pks_ktr, locs_ktr, width_ktr, prom_ktr, heights_ktr] = globalpeaks(metrics.time_series_smoothed_ktr(i,1:min([48+StimulationTimePoint,p.Results.MinLifetime])),5);
    [pks_ktr, locs_ktr, width_ktr, prom_ktr, heights_ktr] = globalpeaks(metrics.time_series_ktr(i,1:min([48+StimulationTimePoint,p.Results.MinLifetime])),10);
    %    [pks_ktr, locs_ktr, width_ktr, prom_ktr, heights_ktr] = globalpeaks(metrics.time_series_ktr(i,1:min([96+StimulationTimePoint,p.Results.MinLifetime])),5);
%    [pks_ktr, locs_ktr, width_ktr, prom_ktr, heights_ktr] = globalpeaks(metrics.time_series_ktr(i,1:min([90,p.Results.MinLifetime])),5);
    % Supress any peaks that are within 4 frames of each other. %Use min difference of 4 for KTR, 6 for NFkB  
    [locs_ktr, order_ktr] = sort(locs_ktr,'ascend');
    pks_ktr = pks_ktr(order_ktr);width_ktr = width_ktr(order_ktr); prom_ktr = prom_ktr(order_ktr); heights_ktr = heights_ktr(order_ktr);
 
    while min(diff(locs_ktr))<4
        tmp = find(diff(locs_ktr)==min(diff(locs_ktr)),1,'first');
        tmp = tmp + (pks_ktr(tmp)>=pks_ktr(tmp+1));
        pks_ktr(tmp) = []; locs_ktr(tmp) = []; width_ktr(tmp) = []; prom_ktr(tmp) = []; heights_ktr(tmp) = [];
    end
    
    pks_ktr(locs_ktr<(StimulationTimePoint + 1)) = [];
    width_ktr(locs_ktr<(StimulationTimePoint + 1)) = [];
    prom_ktr(locs_ktr<(StimulationTimePoint + 1)) = [];
    heights_ktr(locs_ktr<(StimulationTimePoint + 1)) = [];
    locs_ktr(locs_ktr<(StimulationTimePoint + 1)) = [];
%    pks(locs<(StimulationTimePoint + 3)) = []; %this seems to make some early peaks not detected --> remove extra padding
%    locs(locs<(StimulationTimePoint + 3)) = [];

% 20200617 Testing: Filter to remove peaks with amplitudes below 0
%todo test this
    locs_ktr(pks_ktr<= 0) = [];
    width_ktr(pks_ktr<= 0) = [];
    prom_ktr(pks_ktr<= 0) = [];
    heights_ktr(pks_ktr<= 0) = [];
    pks_ktr(pks_ktr<= 0) = []; %pks needs to be filtered after others
    
%20200826 SL Testing: Filter to remove peaks with short peak prominence based on Stdv of baseline
 %testing 3*, 2*
%    locs_ktr(prom_ktr< 2*baseline_stdv_ktr(i)) = [];
 %   width_ktr(prom_ktr< 2*baseline_stdv_ktr(i)) = [];
  %  pks_ktr(prom_ktr< 2*baseline_stdv_ktr(i)) = []; 
   % heights_ktr(prom_ktr< 2*baseline_stdv_ktr(i)) = [];
    %prom_ktr(prom_ktr< 2*baseline_stdv_ktr(i)) = [];%prom pks needs to be filtered after others

    locs_ktr(heights_ktr< 2*baseline_stdv_ktr(i)) = [];
    width_ktr(heights_ktr< 2*baseline_stdv_ktr(i)) = [];
    pks_ktr(heights_ktr< 2*baseline_stdv_ktr(i)) = []; 
    prom_ktr(heights_ktr< 2*baseline_stdv_ktr(i)) = [];%heights pks needs to be filtered after others
    heights_ktr(heights_ktr< 2*baseline_stdv_ktr(i)) = []; %heights pks needs to be filtered after others

%20200923 SL Testing: Filter to remove peaks that are too narrow (based on 'halfheight' width determination in globablpeaks)
    locs_ktr(width_ktr< 2) = [];
    pks_ktr(width_ktr< 2) = []; 
    prom_ktr(width_ktr< 2) = [];
    heights_ktr(width_ktr< 2)= []; 
    width_ktr(width_ktr< 2)= [];%width of pks needs to be filtered last 
    
    
   if ~isempty(locs_ktr)
        metrics.pk1_time_ktr(i) = locs_ktr(1);
        metrics.pk1_amp_ktr(i) = pks_ktr(1);
        metrics.pk1_width_ktr(i) = width_ktr(1);
        metrics.pk1_prom_ktr(i) = prom_ktr(1);
        metrics.pk1_height_ktr(i) = heights_ktr(1);

    end
    if length(locs_ktr)>1
        metrics.pk2_time_ktr(i) = locs_ktr(2);
        metrics.pk2_amp_ktr(i) = pks_ktr(2);
        metrics.pk2_width_ktr(i) = width_ktr(2);
        metrics.pk2_prom_ktr(i) = prom_ktr(2);
        metrics.pk2_height_ktr(i) = heights_ktr(2);

    end
end

metrics.pk1_time_ktr = (metrics.pk1_time_ktr-StimulationTimePoint)/FramesPerHour;
metrics.pk2_time_ktr = (metrics.pk2_time_ktr-StimulationTimePoint)/FramesPerHour;

%% KTR METRICS OF DURATION
% Envelope width: maximum consecutive time above a threshold (envelope must begin within 1st 6 hrs)
smoothed2 = smoothrows(metrics.time_series_ktr,5);
upperThresh = 0.28;
%aux.ktr.thresholds = linspace(0, OnThreshKTR*3, 40);
aux.ktr.thresholds = linspace(0, upperThresh, 25);
metrics.envelope_ktr = zeros(size(metrics.time_series_ktr,1),length(aux.ktr.thresholds));
for j = 1:length(aux.ktr.thresholds)
    thresholded = smoothed2(:,StimulationTimePoint:end)>aux.ktr.thresholds(j); %consider only TP after stimulation
    for i = 1:size(thresholded,1)
        curr = 1;
        idx_start = 1;
        while (curr<size(thresholded,2)) && (idx_start< (6*FramesPerHour))
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
metrics.envelope_ktr = metrics.envelope_ktr/FramesPerHour;


% Number of frames above a given threshold
metrics.duration_ktr = zeros(size(metrics.time_series_ktr,1),length(aux.ktr.thresholds));
%only include TP from stimualtions onwards
for i = 1:length(aux.ktr.thresholds)
    metrics.duration_ktr(:,i) = nansum(smoothed(:,StimulationTimePoint:end)>aux.ktr.thresholds(i),2)/FramesPerHour;
end
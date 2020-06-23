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
addParameter (p, 'GraphLimitsNFkB',[-0.25 8],@isnumeric);
addParameter (p, 'GraphLimitsKTR',[-0.0,0.4],@isnumeric);
expectedFlags = {'on','off'};
addParameter(p,'Verbose','off', @(x) any(validatestring(x,expectedFlags)))
valid_conv = @(x) assert(isnumeric(x)&&(x>=0)&&(length(x)==1),...
    'Convection correction parameter must be single integer >= 0');
addParameter(p,'ConvectionShift',1, valid_conv);
addParameter(p, 'StimulationTimePoint', 13, @isnumeric)
addParameter(p, 'FramesPerHour', 12, @isnumeric)

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
                            'GraphLimitsNFkB', p.Results.GraphLimitsNFkB, 'GraphLimitsKTR', p.Results.GraphLimitsKTR, 'StimulationTimePoint', p.Results.StimulationTimePoint);
   
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

% version including integral below 0  
% 2) integrated activity
% use simple baseline deducted time series
% include time before stimulation (should be close to 0)
%todo check what quantitative effect my new baseline method has, since this now includes values below 0

metrics.integrals_nfkb = cumtrapz(t,metrics.time_series_nfkb,2);

%metrics.integrals_nfkb = nan(size(metrics.time_series_nfkb));
%for i = 1:size(metrics.integrals_nfkb,1)
%    metrics.integrals_nfkb(i,:) = cumtrapz(t,metrics.time_series_nfkb(i,:));
%end

% 3) differentiated activity - use central finite difference
% ? unit: nfkb activity (au)/hour
smoothed = smoothrows(metrics.time_series_nfkb,3);
metrics.derivatives_nfkb = (smoothed(:,3:end) - smoothed(:,1:end-2))/(1/6); %?divided by 1/6 because 2 tp, ie 10 min is 1/6 of an hour?


%% TRIM EVERYBODY to a common length (of "good" sets, current minimum is roughly 21 hrs)
try
    metrics.time_series_nfkb = metrics.time_series_nfkb(:,1:p.Results.TrimFrame);
    metrics.time_series_nfkb_no_base_ded = metrics.time_series_nfkb_no_base_ded(:,1:p.Results.TrimFrame);
    metrics.integrals_nfkb = metrics.integrals_nfkb(:,1:p.Results.TrimFrame);
    metrics.derivatives_nfkb = metrics.derivatives_nfkb(:,1:(p.Results.TrimFrame-2));
    smoothed = smoothed(:,1:p.Results.TrimFrame);
    t = t(1:p.Results.TrimFrame);

catch me
    disp(['Note: vectors too short to cap @ ',num2str(p.Results.TrimFrame),' frames'])
end
%% MISC NFkB METRICS

% 4) Integrals within one-hour windows (0-1, 1-2, 2-3) after stimulation and three hour windows (0-3, 1-4, etc) of activity after stimulation
% use baseline deducted nfkb values
% seems to automatically only include time after stimulation (as long as t is adjusted to have 0 at StimulationTimePoint
max_hr = floor(max(t));
metrics.intwin1_nfkb = nan(size(metrics.time_series_nfkb,1),max_hr);
metrics.intwin3_nfkb = nan(size(metrics.time_series_nfkb,1),max_hr-2);
for i = 1:(max_hr)
    win = t>=(i-1) & t<(i); %eg first window goes from 0 to 1 h post stim (corresponds to TP13 to TP24 with 13 as StimulationTimePoint)
    %double check , shouldnt it go to 25?
    %second window goes from 25 to 36 --> thus act between 24 and 25 is never included?
    metrics.intwin1_nfkb(:,i) = trapz(t(win),metrics.time_series_nfkb(:,win),2);
    if i<= (max_hr-2)
        win = t>=(i-1) & t<(i+2); %eg first window goes from 0 to 3 h post stim (corresponds to TP13 to TP48 with 13 as StimulationTimePoint)
        metrics.intwin3_nfkb(:,i) = trapz(t(win),metrics.time_series_nfkb(:,win),2);
    end
end

%% MAX/MIN metrics
%adjusted to include only time after stimulation
metrics.max_amplitude_nfkb = nanmax(metrics.time_series_nfkb(:,StimulationTimePoint:end),[],2);
metrics.max_integral_nfkb = nanmax(metrics.integrals_nfkb(:,StimulationTimePoint:end),[],2);
metrics.max_derivative_nfkb = nanmax(metrics.derivatives_nfkb(:,StimulationTimePoint:end),[],2);
metrics.min_derivative_nfkb = nanmin(metrics.derivatives_nfkb(:,StimulationTimePoint:end),[],2);

%% ACTIVITY metrics:
%compute responder index using sigma threshold and off times
Wliml = StimulationTimePoint+1; %first/lower time point of window to check for activity
Wlimu = StimulationTimePoint + 48; %last/upper time point of window to check for activity, ie check in the first 4 hours after stimulation
blockLengthThresh = 5; %number of consecutive frames cell needs to pass activity threshold to be considered a responder
baseline_stdv_nfkb = nanstd(metrics.time_series_nfkb(:,1:StimulationTimePoint),0,2);
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
                    %block_length = thresh_stop(1);
                elseif isempty(thresh_stop)
                    %block_length = numel(thresh_idx)+1 - thresh_start(end);
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
%todo write this into another function  
metrics.responders_fraction_nfkb = nnz(metrics.responder_index_nfkb)/numel(metrics.responder_index_nfkb);
%%
% Determination of off_times
% using smoothed values and stdv of non-smoothed basline 
smoothed_by_sigma = smoothed./baseline_stdv_nfkb;
on_array = zeros(size(smoothed(:,Wliml:end)));
%on_array = zeros(size(smoothed));
metrics.off_times_nfkb = zeros(size(smoothed,1),1);
for ii = 1:size(smoothed_by_sigma,1)
    n = 0;
    for jj = Wliml:size(smoothed_by_sigma,2)
  %  for jj = 1: size(smoothed_by_sigma,2)
        if smoothed_by_sigma(ii,jj)> OnThreshNFkB
            n = n+1;
        else
        n = 0;
        end
        %on_array(ii,(jj)) = n;
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
%}
%%
%{
OnThreshNFkB = 1;
window_sz = 14; % ~1 hr windows (on either side of a given timepoint)
thresh = 0.9; % Pct of inactivity allowed in a given window
cutoff_time = 4; % time to look for cell activity before declaring it "off" (hrs)
%Brook's way of determining activity and off times
% ACTIVITY metrics: compute an off-time for all cells
metrics.off_times_nfkb = zeros(size(smoothed,1),1);
inactive = [repmat(nanmin(smoothed(:,1:7),[],2),1,window_sz*2+1),smoothed(:,:),...
    repmat(nanmedian(smoothed(:,(end-window_sz:end)),2),1,window_sz*2)];% creates matrix of min of first couple of TP, repeated 14*2+1 times, all the smoothed rows, and the median of the last 14 rows repeated 14*2 times
% the folliwng somehow 'magically' converts all the numbers that are
% below the activation threshold into a percent inactive, ie 0 is 100% inactive, more active cells are fraction of that
% this happens using a smoothing moving average, but I do not understand how this works?
inactive = smoothrows(inactive<(OnThreshNFkB),(window_sz*2));%smoothes again only those lower than threshold with BIG smoothing window
frontcrop = round(window_sz*2*(1-thresh))+window_sz+1; % do not understand why the threshold appear here
inactive = inactive(:,frontcrop:end); %removes a certain number (18) of the extra columns added above (even though 29 columns were added?)
inactive = inactive(:,1:size(smoothed,2)); %cuts off the end columns added above
inactive(isnan(smoothed)) = nan;

% Find the final time each cell was active
for i = 1:length(metrics.off_times_nfkb)
    active_times = find(inactive(i,:)<thresh);
    if ~isempty(active_times)
        if active_times(1) < (cutoff_time*12) % ignore cells who only turned on after cutoff time hrs.
            metrics.off_times_nfkb(i) = active_times(end);
        end
    end    
end
%metrics.off_times_nfkb = (metrics.off_times_nfkb-1)/12;
metrics.off_times_nfkb = (metrics.off_times_nfkb-StimulationTimePoint)/12;
metrics.off_times_nfkb(metrics.off_times_nfkb<0) = 0;
 metrics.off_times_nfkb_BT_1 = metrics.off_times_nfkb;
%}

%% NFkB OSCILLATION METRICS

% Calculate fourier distribution (via FFT) & power
%todo check if off_pad is needed
off_pad = 12; % Signal time added to trajectory in  FFT calculation (keeps transients from being recorded as osc.)
Fs = 1/300;
depth = max(metrics.off_times_nfkb)*FramesPerHour;
NFFT = 2^nextpow2(depth); % Next power of 2 from chosen depth, can be used to pad input to fft
aux.fft = zeros(size(metrics.time_series_nfkb,1),NFFT/2+1);
aux.freq = Fs/2*linspace(0,1,NFFT/2+1);
aux.power = zeros(size(aux.fft));


for i = 1:size(metrics.time_series_nfkb,1)
    if(metrics.off_times_nfkb(i)>0)
        %adjusted to start with stimulation
        y = metrics.time_series_nfkb(i,StimulationTimePoint:(depth+StimulationTimePoint));
%        y = metrics.time_series_nfkb(i,1:depth);
        off_frame = min([length(y), metrics.off_times_nfkb(i)*FramesPerHour+1+off_pad]); % (Pad w/ 1 extra hr of content)
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
freq_thresh = aux.freq((aux.freq >= (0.35/3600)) & (aux.freq <= (0.7/3600)));
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
  %  [pks, locs] =
  %  globalpeaks(metrics.time_series_nfkb(i,1:min([90,p.Results.MinLifetime])),5);
  %  %20200617 test to include more peaks because of more filtering
    [pks, locs] = globalpeaks(metrics.time_series_nfkb(i,1:min([90,p.Results.MinLifetime])),20);
    % Supress any peaks that are within 6 frames of each other.
    [locs, order] = sort(locs,'ascend');
    pks = pks(order);
    while min(diff(locs))<6
        tmp = find(diff(locs)==min(diff(locs)),1,'first');
        tmp = tmp + (pks(tmp)>=pks(tmp+1));
        pks(tmp) = [];
        locs(tmp) = [];  
    end
    pks(locs<(StimulationTimePoint + 1)) = [];
    locs(locs<(StimulationTimePoint + 1)) = [];
%    pks(locs<(StimulationTimePoint + 3)) = []; %this seems to make some early peaks not detected --> remove extra padding
%    locs(locs<(StimulationTimePoint + 3)) = [];
%    pks(locs<4) = [];
%    locs(locs<4) = [];

% 20200617 Testing to remove peaks with values below 0
%todo test this
    locs(pks<= 0) = [];
    pks(pks<= 0) = []; %pks needs to be filtered after locs
 

%todo include what to do if it is empty, or if pks is empty
   if ~isempty(locs)
        metrics.pk1_time_nfkb(i) = locs(1);
        metrics.pk1_amp_nfkb(i) = pks(1);
    end
    if length(locs)>1
        metrics.pk2_time_nfkb(i) = locs(2);
        metrics.pk2_amp_nfkb(i) = pks(2);
    end
end
metrics.pk1_time_nfkb = (metrics.pk1_time_nfkb-StimulationTimePoint)/FramesPerHour;
metrics.pk2_time_nfkb = (metrics.pk2_time_nfkb-StimulationTimePoint)/FramesPerHour;
%metrics.pk1_time_nfkb = (metrics.pk1_time_nfkb-1)/12;
%metrics.pk2_time_nfkb = (metrics.pk2_time_nfkb-1)/12;


%% NFkB METRICS OF DURATION
% Envelope width: maximum consecutive time above a threshold (envelope must begin within 1st 6 hrs)
smoothed2 = smoothrows(metrics.time_series_nfkb,5);
%aux.thresholds = linspace(0, OnThreshNFkB*3, 40);
upperThresh = 8;
%todo pick proper threshold
%Ade uses only 25 thresholds, not 40
aux.thresholds = linspace(0, upperThresh, 25);
metrics.envelope_nfkb = zeros(size(metrics.time_series_nfkb,1),length(aux.thresholds));
for j = 1:length(aux.thresholds)
    thresholded = smoothed2(:,StimulationTimePoint:end)>aux.thresholds(j);
%    thresholded = smoothed2>aux.thresholds(j);
    for i = 1:size(thresholded,1)
        curr = 1;
        idx_start = 1;
        %todo check if the + stimulationtimepoint is really needed
        while (curr<size(thresholded,2)) && (idx_start< (6*FramesPerHour))
%        while (curr<size(thresholded,2)) && (idx_start< (6*FramesPerHour+StimulationTimePoint))
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
metrics.duration_nfkb = zeros(size(metrics.time_series_nfkb,1),length(aux.thresholds));
for i = 1:length(aux.thresholds)
    metrics.duration_nfkb(:,i) = nansum(smoothed(:,StimulationTimePoint:end)>aux.thresholds(i),2)/FramesPerHour;
%    metrics.duration_nfkb(:,i) = nansum(smoothed>aux.thresholds(i),2)/FramesPerHour;
end
%% NFkB METRICS OF DURATION Using sigma of baseline
%
% Envelope width: maximum consecutive time above a threshold (envelope must begin within 1st 6 hrs after stimulation)
smoothed2 = smoothrows(metrics.time_series_nfkb,5);
%baseline_stdv_nfkb = nanstd(metrics.time_series_nfkb(:,1:StimulationTimePoint),0,2); %already calculated above
smoothed2_by_sigma= smoothed2./baseline_stdv_nfkb;
%todo pick proper threshold
%aux.thresholds = linspace(0, OnThreshNFkB*3, 40);
upperThresh = OnThreshNFkB*12;
aux.thresholds = linspace(0, upperThresh, 25);
metrics.envelope_sigma_nfkb = zeros(size(metrics.time_series_nfkb,1),length(aux.thresholds));
for j = 1:length(aux.thresholds)
    thresholded = smoothed2_by_sigma(:,StimulationTimePoint:end)>aux.thresholds(j);
%    thresholded = smoothed2>aux.thresholds(j);
    for i = 1:size(thresholded,1)
        curr = 1;
        idx_start = 1;

%todo check if the + stimulationtimepoint is really needed
        while (curr<size(thresholded,2)) && (idx_start< (6*FramesPerHour+StimulationTimePoint))
            idx_start = find(thresholded(i,curr:end)==1,1,'first')+curr-1;
            if ~isempty(idx_start)
                idx_stop = find(thresholded(i,idx_start:end)==0,1,'first')+idx_start-1;
                if isempty(idx_stop)
                    idx_stop = find(~isnan(thresholded(i,:)),1,'last');
                end
                if (idx_stop-idx_start) > metrics.envelope_sigma_nfkb(i,j)
                    metrics.envelope_sigma_nfkb(i,j) = (idx_stop-idx_start);
                end
                curr = idx_stop;
            else
                break
            end
        end
    end
end
metrics.envelope_sigma_nfkb = metrics.envelope_sigma_nfkb/FramesPerHour;


% Number of frames above a given threshold
metrics.duration_sigma_nfkb = zeros(size(metrics.time_series_nfkb,1),length(aux.thresholds));
for i = 1:length(aux.thresholds)
    metrics.duration_sigma_nfkb(:,i) = nansum(smoothed_by_sigma(:,StimulationTimePoint:end)>aux.thresholds(i),2)/FramesPerHour;
%    metrics.duration_nfkb(:,i) = nansum(smoothed>aux.thresholds(i),2)/FramesPerHour;
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
%prev integration incl values below 0
    metrics.integrals_ktr = cumtrapz(t,metrics.time_series_ktr,2);
%metrics.integrals_ktr = nan(size(metrics.time_series_ktr));
%for i = 1:size(metrics.integrals_ktr,1)
%    metrics.integrals_ktr(i,:) = cumtrapz(t,metrics.time_series_ktr(i,:));
%end
%}

% 3) differentiated activity - use central finite difference
% use baseline deducted KTR ratios for this, shouldn't make a difference 
smoothed = smoothrows(metrics.time_series_ktr,3);
metrics.derivatives_ktr = (smoothed(:,3:end) - smoothed(:,1:end-2))/(1/6);


%% TRIM EVERYBODY to a common length (of "good" sets, current minimum is roughly 21 hrs)
try
    metrics.time_series_ktr = metrics.time_series_ktr(:,1:p.Results.TrimFrame);
    metrics.time_series_ktr_no_base_ded = metrics.time_series_ktr_no_base_ded(:,1:p.Results.TrimFrame);
    metrics.integrals_ktr = metrics.integrals_ktr(:,1:p.Results.TrimFrame);
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
metrics.max_integral_ktr = nanmax(metrics.integrals_ktr(:,StimulationTimePoint:end),[],2);
metrics.max_derivative_ktr = nanmax(metrics.derivatives_ktr(:,StimulationTimePoint:end),[],2);
metrics.min_derivative_ktr = nanmin(metrics.derivatives_ktr(:,StimulationTimePoint:end),[],2);

%% ACTIVITY metrics: compute an off-time for all cells

%compute responder index using sigma threshold and off times
Wliml = StimulationTimePoint+1; %first/lower time point of window to check for activity
Wlimu = StimulationTimePoint + 36; %last/upper time point of window to check for activity, ie check in the first 4 hours after stimulation
blockLengthThresh = 3; %number of consecutive frames cell needs to pass activity threshold to be considered a responder
baseline_stdv_ktr = nanstd(metrics.time_series_ktr(:,1:StimulationTimePoint),0,2);
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
                    %block_length = thresh_stop(1);
                elseif isempty(thresh_stop)
                    %block_length = numel(thresh_idx)+1 - thresh_start(end);
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
%todo write this into another function      
metrics.responders_fraction_ktr = nnz(metrics.responder_index_ktr)/numel(metrics.responder_index_ktr);
%%
%
%Attempt at simplistic determination of off_times
%using smoothed data and stdv of non-smoothedbasline
smoothed_by_sigma = smoothed./baseline_stdv_ktr;
on_array = zeros(size(smoothed(:,Wliml:end)));
% on_array = zeros(size(smoothed));
metrics.off_times_ktr = zeros(size(smoothed,1),1);
for ii = 1:size(smoothed_by_sigma,1)
    n = 0;
    for jj = Wliml:size(smoothed_by_sigma,2)
%    for jj = 1:size(smoothed_by_sigma,2)
        if smoothed_by_sigma(ii,jj)> OnThreshKTR
            n = n+1;
        else
        n = 0;
        end
        %on_array(ii,jj) = n;
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
%metrics.off_times_ktr = (metrics.off_times_ktr-StimulationTimePoint)/FramesPerHour;
metrics.off_times_ktr(metrics.off_times_ktr<0) = 0;
%}
%%
%{
OnThreshKTR = 0.1;
window_sz = 14; % ~1 hr windows (on either side of a given timepoint)
thresh = 0.9; % Pct of inactivity allowed in a given window
cutoff_time = 4; % time to look for cell activity before declaring it "off" (hrs)
%Brook's way of determining activity and off times
% ACTIVITY metrics: compute an off-time for all cells
metrics.off_times_ktr = zeros(size(smoothed,1),1);
inactive = [repmat(nanmin(smoothed(:,1:7),[],2),1,window_sz*2+1),smoothed(:,:),...
    repmat(nanmedian(smoothed(:,(end-window_sz:end)),2),1,window_sz*2)];% creates matrix of min of first couple of TP, repeated 14*2+1 times, all the smoothed rows, and the median of the last 14 rows repeated 14*2 times
% the folliwng somehow 'magically' converts all the numbers that are
% below the activation threshold into a percent inactive, ie 0 is 100% inactive, more active cells are fraction of that
% this happens using a smoothing moving average, but I do not understand how this works?
inactive = smoothrows(inactive<(OnThreshKTR),(window_sz*2));%smoothes again only those lower than threshold with BIG smoothing window
frontcrop = round(window_sz*2*(1-thresh))+window_sz+1; % do not understand why the threshold appear here
inactive = inactive(:,frontcrop:end); %removes a certain number (18) of the extra columns added above (even though 29 columns were added?)
inactive = inactive(:,1:size(smoothed,2)); %cuts off the end columns added above
inactive(isnan(smoothed)) = nan;

% Find the final time each cell was active
for i = 1:length(metrics.off_times_ktr)
    active_times = find(inactive(i,:)<thresh);
    if ~isempty(active_times)
        if active_times(1) < (cutoff_time*FramesPerHour) % ignore cells who only turned on after cutoff time hrs.
            metrics.off_times_ktr(i) = active_times(end);
        end
    end    
end
%metrics.off_times_ktr = (metrics.off_times_ktr-1)/FramesPerHour;
metrics.off_times_ktr = (metrics.off_times_ktr-StimulationTimePoint)/FramesPerHour;
metrics.off_times_ktr(metrics.off_times_ktr<0) = 0;
 metrics.off_times_ktr_BT_1 = metrics.off_times_ktr;
%}


%% KTR OSCILLATION METRICS
% Calculate fourier distribution (via FFT) & power
%todo check of off_pad is needed with my new offtimes, etc
off_pad = FramesPerHour; % Signal time added to trajectory in  FFT calculation (keeps transients from being recorded as osc.)
Fs = 1/300;
depth = max(metrics.off_times_ktr)*FramesPerHour;
NFFT = 2^nextpow2(depth); % Next power of 2 from chosen depth
aux.fft = zeros(size(metrics.time_series_ktr,1),NFFT/2+1);
aux.freq = Fs/2*linspace(0,1,NFFT/2+1);
aux.power = zeros(size(aux.fft));

for i = 1:size(metrics.time_series_ktr,1)
    if(metrics.off_times_ktr(i)>0) %only calculate peakfreq if offtime is not 0
        y = metrics.time_series_ktr(i,StimulationTimePoint:(depth+StimulationTimePoint));
%       y = metrics.time_series_ktr(i,1:depth);
        off_frame = min([length(y), metrics.off_times_ktr(i)*FramesPerHour+1+off_pad]); % (Pad w/ 1 extra hr of content)
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
%todo check if these thresholds need adjusting for KTR
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
%}

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
    pks(locs<(StimulationTimePoint + 1)) = [];
    locs(locs<(StimulationTimePoint + 1)) = [];
%    pks(locs<(StimulationTimePoint + 3)) = [];
%    locs(locs<(StimulationTimePoint + 3)) = [];
%    pks(locs<4) = [];
%    locs(locs<4) = [];

% 20200617 Testing to remove peaks with values below 0
%todo test this
%todo inlcude this in KTR section too
    locs(pks<= 0) = [];
    pks(pks<= 0) = []; %pks needs to be filtered after locs
    
    if ~isempty(locs)
        metrics.pk1_time_ktr(i) = locs(1);
        metrics.pk1_amp_ktr(i) = pks(1);
    end
    if length(locs)>1
        metrics.pk2_time_ktr(i) = locs(2);
        metrics.pk2_amp_ktr(i) = pks(2);
    end
end
metrics.pk1_time_ktr = (metrics.pk1_time_ktr-StimulationTimePoint)/FramesPerHour;
metrics.pk2_time_ktr = (metrics.pk2_time_ktr-StimulationTimePoint)/FramesPerHour;
%metrics.pk1_time_ktr = (metrics.pk1_time_ktr-1)/FramesPerHour;
%metrics.pk2_time_ktr = (metrics.pk2_time_ktr-1)/FramesPerHour;


%% KTR METRICS OF DURATION
% Envelope width: maximum consecutive time above a threshold (envelope must begin within 1st 6 hrs)
smoothed2 = smoothrows(metrics.time_series_ktr,5);
%fixed use of OnThresh, but using other threshold (0.8 is lower than most max amplitude across a couple of exp
%todo pick proper threshold
%aux.thresholds = linspace(0, OnThreshKTR*3, 40);
upperThresh = 0.3;
%20200606 temp make threshold lower to test CC calc (make sure less 0) in metrics
%upperThresh = 0.4;
aux.thresholds = linspace(0, upperThresh, 25);
metrics.envelope_ktr = zeros(size(metrics.time_series_ktr,1),length(aux.thresholds));
for j = 1:length(aux.thresholds)
    thresholded = smoothed2(:,StimulationTimePoint:end)>aux.thresholds(j); %consider only TP after stimulation
    for i = 1:size(thresholded,1)
        curr = 1;
        idx_start = 1;
%todo check if the + stimulationtimepoint is really needed
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
metrics.duration_ktr = zeros(size(metrics.time_series_ktr,1),length(aux.thresholds));
%only include TP from stimualtions onwards
for i = 1:length(aux.thresholds)
    metrics.duration_ktr(:,i) = nansum(smoothed(:,StimulationTimePoint:end)>aux.thresholds(i),2)/FramesPerHour;
end

smoothed2_by_sigma= smoothed2./baseline_stdv_nfkb;
%% KTR METRICS OF DURATION using sigma thresholds
% Envelope width: maximum consecutive time above a threshold (envelope must begin within 1st 6 hrs)
smoothed2 = smoothrows(metrics.time_series_ktr,5);
smoothed2_by_sigma= smoothed2./baseline_stdv_ktr;
%fixed use of OnThresh, but using other threshold (0.8 is lower than most max amplitude across a couple of exp
%todo pick proper threshold
%aux.thresholds = linspace(0, OnThreshKTR*3, 40);
%SL 20200606 temp change in threshold to test CC calc
%upperThresh = OnThreshKTR*5;
upperThresh = OnThreshKTR*4;
aux.thresholds = linspace(0, upperThresh, 25);
metrics.envelope_sigma_ktr = zeros(size(metrics.time_series_ktr,1),length(aux.thresholds));
for j = 1:length(aux.thresholds)
    thresholded = smoothed2_by_sigma(:,StimulationTimePoint:end)>aux.thresholds(j); %consider only TP after stimulation
    for i = 1:size(thresholded,1)
        curr = 1;
        idx_start = 1;
 %todo check if the + stimulationtimepoint is really needed
        while (curr<size(thresholded,2)) && (idx_start< (6*FramesPerHour))
            idx_start = find(thresholded(i,curr:end)==1,1,'first')+curr-1;
            if ~isempty(idx_start)
                idx_stop = find(thresholded(i,idx_start:end)==0,1,'first')+idx_start-1;
                if isempty(idx_stop)
                    idx_stop = find(~isnan(thresholded(i,:)),1,'last');
                end
                if (idx_stop-idx_start) > metrics.envelope_sigma_ktr(i,j)
                    metrics.envelope_sigma_ktr(i,j) = (idx_stop-idx_start);
                end
                curr = idx_stop;
            else
                break
            end
        end
    end
end
metrics.envelope_sigma_ktr = metrics.envelope_sigma_ktr/FramesPerHour;


% Number of frames above a given threshold
metrics.duration_sigma_ktr = zeros(size(metrics.time_series_ktr,1),length(aux.thresholds));
%only include TP from stimualtions onwards
for i = 1:length(aux.thresholds)
    metrics.duration_sigma_ktr(:,i) = nansum(smoothed_by_sigma(:,StimulationTimePoint:end)>aux.thresholds(i),2)/FramesPerHour;
end

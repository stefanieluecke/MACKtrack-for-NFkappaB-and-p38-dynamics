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
addParameter(p, 'MinSize', 90, @isnumeric); 
addParameter(p,'MinLifetime',109, @isnumeric);
addParameter(p,'TrimFrame',157, @isnumeric);
addParameter (p, 'GraphLimitsNFkB',[-0.25 7],@isnumeric);
addParameter (p, 'GraphLimitsKTR',[-0.02,0.35],@isnumeric);
expectedFlags = {'on','off'};
addParameter(p,'Verbose','off', @(x) any(validatestring(x,expectedFlags)))
addParameter(p, 'StimulationTimePoint', 13, @isnumeric)
addParameter(p, 'FramesPerHour', 12, @isnumeric)
addParameter(p, 'NFkBBackgroundAdjustment', 'on',@(x) any(validatestring(x,expectedFlags))) %option to turn off NFkB fluorescence distribution adjustment
addParameter(p,'NFkBBaselineDeduction', 'on', @(x) any(validatestring(x,expectedFlags))) %option to turn off NFkB baseline deduction
addParameter(p,'NFkBBaselineAdjustment', 'on', @(x) any(validatestring(x,expectedFlags))) %option to turn off adjusment of NFkB trajectories with correction factor for fluorescence drop derived from Mock experiments
addParameter(p, 'BrooksBaseline', 'off', @(x) any(validatestring(x,expectedFlags)))
addParameter(p,'KTRBaselineDeduction', 'on', @(x) any(validatestring(x,expectedFlags))) %option to turn off NFkB baseline deduction
addParameter(p, 'KTRBaselineAdjustment', 'on',@(x) any(validatestring(x,expectedFlags))) %option to turn off KTR fluorescence distribution adjustment
addParameter(p, 'IncludeKTR', 'on',@(x)any(validatestring(x, expectedFlags)));

parse(p,id, varargin{:})

%% INITIALIZATION. Load and process data

OnThreshNFkB = p.Results.OnThreshNFkB; % Minimum activity required for cell to register as 'on'
OnThreshKTR = p.Results.OnThreshKTR; % Minimum activity required for cell to register as 'on'
StimulationTimePoint = p.Results.StimulationTimePoint;
MinSize = p.Results.MinSize; 
MinLifetime = p.Results.MinLifetime; 
FramesPerHour = p.Results.FramesPerHour;

[graph, info, measure] = filter_nfkb_ktr_ratio(id,'MinLifetime',MinLifetime,...
                             'OnThreshNFkB',OnThreshNFkB,...
                            'OnThreshKTR',OnThreshKTR,'MinSize', MinSize, 'Verbose', p.Results.Verbose,...
                            'GraphLimitsNFkB', p.Results.GraphLimitsNFkB, 'GraphLimitsKTR', p.Results.GraphLimitsKTR, 'StimulationTimePoint', p.Results.StimulationTimePoint,...
                            'FramesPerHour', p.Results.FramesPerHour, 'NFkBBaselineDeduction', p.Results.NFkBBaselineDeduction, 'NFkBBackgroundAdjustment',p.Results.NFkBBackgroundAdjustment,'NFkBBaselineAdjustment', p.Results.NFkBBaselineAdjustment,... 
                            'BrooksBaseline', p.Results.BrooksBaseline, 'KTRBaselineDeduction', p.Results.KTRBaselineDeduction,'KTRBaselineAdjustment', p.Results.KTRBaselineAdjustment, 'IncludeKTR', p.Results.IncludeKTR);
   
graph.var_nfkb = graph.var_nfkb(:,1:min(p.Results.TrimFrame, size(graph.var_nfkb,2))); 
graph.var_nfkb_no_base_ded = graph.var_nfkb_no_base_ded(:,1:min(p.Results.TrimFrame, size(graph.var_nfkb_no_base_ded,2)));

if strcmpi(p.Results.IncludeKTR,'on')
    graph.var_ktr = graph.var_ktr(:,1:min(p.Results.TrimFrame, size(graph.var_ktr,2)));
    graph.var_ktr_no_base_ded = graph.var_ktr_no_base_ded(:,1:min(p.Results.TrimFrame, size(graph.var_ktr_no_base_ded,2)));
end
graph.t = graph.t(1:size(graph.var_nfkb,2));

graph.opt_nfkb = maketicks(graph.t,info.GraphLimitsNFkB,0); 
graph.opt_nfkb.Name = 'NF\kappaB activation'; 

if strcmpi(p.Results.IncludeKTR,'on')
    graph.opt_ktr = maketicks(graph.t,info.GraphLimitsKTR,0); 
    graph.opt_ktr.Name = 'Kinase activation'; 
end


if ~ismember ('OnThreshNFkB',p.UsingDefaults)
    OnThreshNFkB = info.OnThreshNFkB;
end
if strcmpi(p.Results.IncludeKTR,'on')
    if ~ismember ('OnThreshKTR',p.UsingDefaults)
        OnThreshKTR = info.OnThreshKTR;
    end
end
%% NFkB METRICS
%% BASIC NFkB METRICS: TIME SERIES, BASELINE, DERIVATIVE, INTEGRAL
% Basic time series. Interpolate over "normal" interval (12 frames per hr) if required
    % use baseline deducted NFkB, also generate non-baseline deducted NFkB here
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

% Include baseline calculated in filter function into metrics
    metrics.baseline_nfkb = info.nfkb_baseline;

% integrated activity
    % includes below 0 activity 
    % start with stimulation time point
    metrics.integrals_nfkb = cumtrapz(t(StimulationTimePoint:end),metrics.time_series_nfkb(:,StimulationTimePoint:end),2);

% total activity = 'maximum integral'
    metrics.max_integral_nfkb = nanmax(metrics.integrals_nfkb,[],2);

% integrals above 0 only
    integrals_pos = nan(size(metrics.integrals_nfkb ));
    integrals_pos(:,1) = metrics.integrals_nfkb (:,1);
    for i = 1:size(integrals_pos,1)
        for j = 2:size(metrics.integrals_nfkb ,2)
            if metrics.integrals_nfkb (i,j)>= metrics.integrals_nfkb (i,j-1)
             integrals_pos(i,j) =  integrals_pos(i,j-1) + (metrics.integrals_nfkb (i,j) - metrics.integrals_nfkb (i,j-1));
            elseif metrics.integrals_nfkb (i,j)< metrics.integrals_nfkb (i,j-1)
             integrals_pos(i,j) =  integrals_pos(i,j-1);
            end
        end
    end

    metrics.integrals_pos_nfkb =integrals_pos; 
% total activity above 0= 'maximum integral'
    metrics.max_pos_integral_nfkb = nanmax(integrals_pos,[],2);
% time2XMaxActivity (within 1st 8 hours or less)
    endFrame = min(96, size(metrics.integrals_nfkb,2));
    P25MaxIntegral = nanmax(integrals_pos(:,1:endFrame),[],2)*0.25; % half max integral
    P50MaxIntegral = nanmax(integrals_pos(:,1:endFrame),[],2)*0.5; % half max integral
    P75MaxIntegral = nanmax(integrals_pos(:,1:endFrame),[],2)*0.75; % half max integral

    distances25 = abs(integrals_pos(:,1:endFrame)- P25MaxIntegral);
    distances50 = abs(integrals_pos(:,1:endFrame)- P50MaxIntegral);
    distances75 = abs(integrals_pos(:,1:endFrame)- P75MaxIntegral);
    [~, idx25] = nanmin(distances25,[],2);
    [~, idx50] = nanmin(distances50,[],2);
    [~, idx75] = nanmin(distances75,[],2);
    idx25(idx25==1) = NaN;
    idx50(idx50==1) = NaN;
    idx75(idx75==1) = NaN;
    metrics.time2XmaxPosInt_nfkb = nan(size(metrics.integrals_pos_nfkb, 1), 3);
    metrics.time2XmaxPosInt_nfkb(:,1)= (idx25-1)/FramesPerHour;
    metrics.time2XmaxPosInt_nfkb(:,2)= (idx50-1)/FramesPerHour;
    metrics.time2XmaxPosInt_nfkb(:,3)= (idx75-1)/FramesPerHour;

%%
% Differentiated activity - use central finite difference
    % use StimulationTimePoint as start
    % end after 4h, i.e. 61 TPs (12 baseline, 1 start, 48 2h stim)
    smoothed = smoothdata(fillmissing(metrics.time_series_nfkb,'linear', 2,'EndValues','extrap'),2,'sgolay', 6); 
%
    metrics.derivatives_nfkb = (smoothed(:,StimulationTimePoint+2:62) - smoothed(:,StimulationTimePoint:60))/(1/6); %divided by 1/6 because 2 tp, ie 10 min is 1/6 of an hour?
%}
%% TRIM EVERYBODY to a common length
try
    metrics.time_series_nfkb = metrics.time_series_nfkb(:,1:p.Results.TrimFrame);
    metrics.time_series_nfkb_no_base_ded = metrics.time_series_nfkb_no_base_ded(:,1:p.Results.TrimFrame);
    metrics.time_series_ps_nfkb = metrics.time_series_nfkb(:,StimulationTimePoint:p.Results.TrimFrame);
    metrics.integrals_nfkb = metrics.integrals_nfkb(:,1:(p.Results.TrimFrame-StimulationTimePoint));
    metrics.integrals_pos_nfkb = metrics.integrals_pos_nfkb(:,1:(p.Results.TrimFrame-StimulationTimePoint));
%    metrics.derivatives_nfkb = metrics.derivatives_nfkb(:,1:(p.Results.TrimFrame-2-StimulationTimePoint));
    smoothed = smoothed(:,1:p.Results.TrimFrame);
    t = t(1:p.Results.TrimFrame);
catch me
    disp(['Note: vectors too short to cap @ ',num2str(p.Results.TrimFrame),' frames'])
end
    metrics.time_series_smoothed_nfkb = smoothed;
%% 

% Integrals within 30-min windows (0-0.5, 0.5-1, 1-1.5), one-hour windows (0-1, 1-2, 2-3) after stimulation and three hour windows (0-3, 1-4, etc) of activity after stimulation
    % automatically only includes time after stimulation (as long as t is adjusted to have 0 at StimulationTimePoint
    max_hr = floor(max(t));
    metrics.intwin1_nfkb = nan(size(metrics.time_series_nfkb,1),max_hr);
    metrics.intwin3_nfkb = nan(size(metrics.time_series_nfkb,1),max_hr-2);
    for i = 1:(max_hr)
        win = t>=(i-1) & t<=(i); %eg first window goes from 0 to 1 h post stim
        metrics.intwin1_nfkb(:,i) = trapz(t(win),metrics.time_series_nfkb(:,win),2);
        if i<= (max_hr-2)
            win = t>=(i-1) & t<=(i+2); %eg first window goes from 0 to 3 h post stim 
            metrics.intwin3_nfkb(:,i) = trapz(t(win),metrics.time_series_nfkb(:,win),2);
        end
    end
    max_half_hr = floor(max(t)*2)/2;
    metrics.intwin_p5_nfkb = nan(size(metrics.time_series_nfkb,1),max_half_hr/0.5);
    counter = 1;
    for i = 0.5:0.5:max_half_hr
        win = t>=(i-0.5) & t<=(i); %eg first window goes from 0 to 0.5 h post stim
        metrics.intwin_p5_nfkb(:,counter) = trapz(t(win),metrics.time_series_nfkb(:,win),2);
        counter = counter + 1;
    end

% Difference between hour 0-1 integral and hour 1-2 integral
%
    metrics.phase_diff1_nfkb = metrics.intwin1_nfkb(:, 1)-metrics.intwin1_nfkb(:, 2);
    metrics.phase_ratio1_nfkb = metrics.intwin1_nfkb(:, 1)./metrics.intwin1_nfkb(:, 2);

% Difference between hour 0-3 integral and hour 3-6 integral
    metrics.phase_diff3_nfkb= metrics.intwin3_nfkb(:, 1)-metrics.intwin3_nfkb(:, 4);
    metrics.phase_ratio3_nfkb= metrics.intwin3_nfkb(:, 1)./metrics.intwin3_nfkb(:, 4);
%}
%% AMPLITUDE metrics
    %include only time after stimulation
    %include only first 4 h
    if size(metrics.time_series_nfkb,2)>=48+StimulationTimePoint
        metrics.max_amp_nfkb = nanmax(metrics.time_series_nfkb(:,StimulationTimePoint:48+StimulationTimePoint),[],2);
    end

    metrics.max2minDiff_nfkb = peak2peak(smoothed(:,StimulationTimePoint:end),2);
%% SIGNAL STATS metrics
    metrics.peak2rms_nfkb     =peak2rms(smoothed(:,StimulationTimePoint:end),2);     
    metrics.rms_nfkb          =rms(smoothed(:,StimulationTimePoint:end),2);
    metrics.mean_movmad_nfkb  =mean( movmad(smoothed(:,StimulationTimePoint:end),3,2),2);
    metrics.mean_movvar_nfkb       =mean( movvar(smoothed(:,StimulationTimePoint:end),3,[],2),2);

%% ACTIVITY metrics:
% Responder status using sigma threshold and off times
% activity has to start within 4 hours
    Wliml = StimulationTimePoint+1; %first/lower time point of window to check for activity
    if size(metrics.time_series_nfkb,2)>48+StimulationTimePoint
    
        Wlimu = StimulationTimePoint + 48; %last/upper time point of window to check for activity, ie check in the first 4 hours after stimulation
    else
        Wlimu = size(metrics.time_series_nfkb,2);
    end
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
                metrics.responder_status_nfkb(jj,1) = any(block_length>=blockLengthThresh, 2);
                block_length_readout(jj, 1:numel(block_length)) = block_length;
    end
    
    metrics.responders_fraction_nfkb = nnz(metrics.responder_status_nfkb)/numel(metrics.responder_status_nfkb);
%%
% Off_times
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
    metrics.off_times_nfkb = metrics.off_times_nfkb/FramesPerHour; %adjustment for stimulation time point not necessary here, because we only start counting frames after stimulation (winl)
    metrics.off_times_nfkb(metrics.off_times_nfkb<0) = 0;

%% NFkB METRICS OF DURATION
    % Envelope width: maximum consecutive time above a threshold (envelope must begin within 1st 6 hrs)
    upperThresh = 7.32;
    aux.nfkb.thresholds = linspace(0, upperThresh, 5);
    metrics.envelope_nfkb = zeros(size(metrics.time_series_nfkb,1),length(aux.nfkb.thresholds));
    for j = 1:length(aux.nfkb.thresholds)
        thresholded = smoothed(:,StimulationTimePoint:end)>aux.nfkb.thresholds(j);
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
%% NFkB OSCILLATION METRICS
%
% Calculate fourier distribution (via FFT) & power
    off_pad = 12; % Signal time added to trajectory in  FFT calculation (keeps transients from being recorded as osc.)
    Fs = 1/300;
    depth = max(metrics.off_times_nfkb)*FramesPerHour;
    NFFT = 2^nextpow2(depth); % Next power of 2 from chosen depth, can be used to pad input to fft
    aux.nfkb.fft = zeros(size(metrics.time_series_nfkb,1),NFFT/2+1);
    aux.nfkb.freq = Fs/2*linspace(0,1,NFFT/2+1);
    aux.nfkb.power = zeros(size(aux.nfkb.fft));
    
    for i = 1:size(metrics.time_series_nfkb,1)
        if(metrics.off_times_nfkb(i)>0) %only calculate peakfreq if offtime is not 0
            %adjusted to start with stimulation
            y = metrics.time_series_nfkb(i,StimulationTimePoint:(depth+StimulationTimePoint));
            off_frame = min([length(y), metrics.off_times_nfkb(i)*FramesPerHour+1+off_pad]); % (Pad w/ 1 extra hr of content)
            y(off_frame:end) = nan;
            y(isnan(y)) = [];
            y = y-nanmean(y);
            if ~isempty(y)
                Y = fft(y,NFFT);
                aux.nfkb.fft(i,:) = abs(Y(1:NFFT/2+1));
                 % Normalization (Plancherel's theorem) to make power peaks comaparable
                aux.nfkb.power(i,:) = abs(Y(1:NFFT/2+1).^2)./sum(abs(y).^2);
    
            end
        end
    end

% Original Peakfreq metric using aux.power
%{
% Find the point of peak (secondary) power
metrics.peakfreq_nfkb = nan(size(aux.nfkb.power,1),1);
for i =1:size(metrics.time_series_nfkb,1)
    [pks,locs] = globalpeaks(aux.nfkb.power(i,:),2);
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
%}
% New Peakfreq metric using normalized aux.power
    % Find the point of peak (secondary) power
    metrics.peakfreq_nfkb = nan(size(aux.nfkb.power,1),1);
    for i =1:size(metrics.time_series_nfkb,1)
        [pks,locs] = globalpeaks(aux.nfkb.power(i,:),2);
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
            metrics.peakfreq_nfkb(i) = 3600*aux.nfkb.freq(idx);
        elseif ~isempty(locs)
             metrics.peakfreq_nfkb(i) = 3600*aux.nfkb.freq(max([locs,3]));
        else
            metrics.peakfreq_nfkb(i) = 3600*aux.nfkb.freq(1);
        end
    end

%%
% Total oscillatory content of particular cells (using thresholds from 0.35 to 0.7 hrs^(-1))
    freq_thresh = aux.nfkb.freq((aux.nfkb.freq >= (0.35/3600)) & (aux.nfkb.freq <= (0.7/3600)));
%{
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
%}

    %using normalized power
    metrics.oscfrac_nfkb = nan(size(aux.nfkb.power,1),length(freq_thresh));
    for j = 1:length(freq_thresh)
        for i =1:size(metrics.time_series_nfkb,1)
            metrics.oscfrac_nfkb(i,j) = nansum(aux.nfkb.power(i,((aux.nfkb.freq >= freq_thresh(j))&(aux.nfkb.power(i,:)>= 6.5)))) /nansum(aux.nfkb.power(i,aux.nfkb.power(i,:)>= 6.5));
    %        metrics.oscfrac_nfkb(i,j) = nansum(aux.nfkb.power(i,aux.nfkb.freq >= freq_thresh(j))) /nansum(aux.nfkb.power(i,:));
            if isnan(metrics.oscfrac_nfkb(i,j))
                metrics.oscfrac_nfkb(i,j) = 0;
            end
        end
    end
 
% NFkB Oscillation power 
    Fs=12; 
    freq_range =[0.35 1];
    n = size(smoothed(:,StimulationTimePoint:end),2); 
    [pwr,fq]=pwelch(smoothed(:,StimulationTimePoint:end)',n,10,256,Fs,'one-sided','psd'); %window size is full length of data? with 10 overlap, 256 is nfft
    fq_nfkb       =fq';
    %normalize power/psd to 1
    psd_nfkb      =transpose(pwr./sum(pwr,1));        
    %oscpower, also referred to bandpower
    pwr = transpose(psd_nfkb) ; fq = transpose(fq_nfkb);% 
    %todo see if freq_range etc need any adjustments esp for KTR
    bp= bandpower(pwr,fq,freq_range, 'psd')';
    metrics.oscpower_nfkb =bp;

% 99% oscillation bandwidth 
    metrics.oscbandwidth_nfkb     =obw(smoothed(:,StimulationTimePoint:end)',Fs)';

%}
%% NFkB METRICS OF AMPLITUDE AND TIMING
%
% 1st peak time/amplitude/prominence/width
    % look for 1st peak in 4h after stimulation
    pk_feats = {'pk1_amp_nfkb', 'pk1_time_nfkb', 'pk1_width_nfkb', 'pk1_prom_nfkb'};
    for i=1:length(pk_feats)
        metrics.(pk_feats{i}) = nan(size(metrics.time_series_nfkb,1),1);
    end
    
    for i = 1:size(metrics.pk1_time_nfkb,1)
        [pks_nfkb, locs_nfkb, width_nfkb, prom_nfkb] = findpeaks(smoothed(i,StimulationTimePoint:(48+StimulationTimePoint)),'NPeaks',8, 'SortStr', 'none', 'WidthReference', 'halfheight', 'MinPeakWidth', 2, 'MinPeakHeight', 0,'MinPeakProminence',1*baseline_stdv_nfkb(i),'Annotate', 'extent');
    %   findpeaks(smoothed(i,StimulationTimePoint:(48+StimulationTimePoint)),'NPeaks',8, 'SortStr', 'none', 'WidthReference', 'halfheight','MinPeakWidth', 2, 'MinPeakHeight', 0,'MinPeakProminence',1*baseline_stdv_nfkb(i),'Annotate', 'extent');
       if ~isempty(locs_nfkb)
            metrics.pk1_time_nfkb(i) = locs_nfkb(1);
            metrics.pk1_amp_nfkb(i) = pks_nfkb(1);
            metrics.pk1_width_nfkb(i) = width_nfkb(1);
            metrics.pk1_prom_nfkb(i) = prom_nfkb(1);
        end
    end
    metrics.pk1_time_nfkb = ((metrics.pk1_time_nfkb)-1)/FramesPerHour;
%% 
% NFkB METRICS OF PERC PK1 AMP DURATION
% Envelope width: maximum consecutive time above a threshold (envelope must begin within 1st 6 hrs)
    % threshold based on percent pk1 amp for each cell
    upperThresh = 0.8;
    aux.nfkb.thresholds = linspace(0.2, upperThresh, 5);
    metrics.envelope_perc_pk1amp_nfkb = zeros(size(metrics.time_series_nfkb,1),length(aux.nfkb.thresholds));
    for j = 1:length(aux.nfkb.thresholds)
        thresholded = smoothed(:,StimulationTimePoint:end)>aux.nfkb.thresholds(j).*metrics.pk1_amp_nfkb;
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
                    if (idx_stop-idx_start) > metrics.envelope_perc_pk1amp_nfkb(i,j)
                        metrics.envelope_perc_pk1amp_nfkb(i,j) = (idx_stop-idx_start);
                        metrics.starttime_env_perc_pk1amp_nfkb(i,j) = idx_start;
    
                    end
                    curr = idx_stop;
                else
                    break
                end
    
            end
            if metrics.envelope_perc_pk1amp_nfkb(i,j) == 0
            metrics.starttime_env_perc_pk1amp_nfkb(i,j) = NaN;
            end
        end
    end
    metrics.envelope_perc_pk1amp_nfkb = metrics.envelope_perc_pk1amp_nfkb/FramesPerHour;
    metrics.starttime_env_perc_pk1amp_nfkb = metrics.starttime_env_perc_pk1amp_nfkb/FramesPerHour;

% Number of frames above a given threshold
    metrics.duration_perc_pk1amp_nfkb = zeros(size(metrics.time_series_nfkb,1),length(aux.nfkb.thresholds));
    for i = 1:length(aux.nfkb.thresholds)
        metrics.duration_perc_pk1amp_nfkb(:,i) = nansum(smoothed(:,StimulationTimePoint:end)>aux.nfkb.thresholds(i).*metrics.pk1_amp_nfkb,2)/FramesPerHour;
    end

%% NFkB Speed metrics
% Max derivative before pk1, min derivative within 1 h after pk1, up and down time to half max pk1 amp
    pk1_frame = metrics.pk1_time_nfkb *FramesPerHour;
    metrics.max_derivative_pk1_nfkb = nan(size(metrics.pk1_time_nfkb));
    max_derivative_pk1_frame = nan(size(metrics.pk1_time_nfkb));
    metrics.min_derivative_pk1_nfkb = nan(size(metrics.pk1_time_nfkb));
    min_derivative_pk1_frame = nan(size(metrics.pk1_time_nfkb));
    %{
    for i = 1:numel(metrics.max_derivative_pk1_nfkb )
        if ~isnan(metrics.pk1_time_nfkb(i)) &&  pk1_frame(i) <= 36
                [metrics.max_derivative_pk1_nfkb(i), max_derivative_pk1_frame(i)] = nanmax(metrics.derivatives_nfkb(i,1:pk1_frame(i)),[],2);
                [metrics.min_derivative_pk1_nfkb(i), min_derivative_pk1_frame(i)] = nanmin(metrics.derivatives_nfkb(i,pk1_frame(i):pk1_frame(i)+FramesPerHour),[],2);
        end
    end
    max_pk1_speed_time_nfkb = max_derivative_pk1_frame./FramesPerHour;
    min_pk1_speed_time_nfkb = min_derivative_pk1_frame./FramesPerHour;
    %}
    metrics.timeUp2halfAmp_nfkb= nan(size(metrics.pk1_time_nfkb));
    metrics.timeDown2halfAmp_nfkb= nan(size(metrics.pk1_time_nfkb));
    half_amp = metrics.pk1_amp_nfkb./2;
    cross_halfAmp = smoothed(:,StimulationTimePoint:end)-half_amp>=0;
    [~, timeUp2halfAmp] = nanmax(cross_halfAmp, [], 2);
    timeDown2halfAmp= nan(size(pk1_frame));
    for idx=1:length(metrics.pk1_time_nfkb)
        if ~isnan(metrics.pk1_time_nfkb(idx))
            cross_halfAmp = smoothed(idx, StimulationTimePoint+pk1_frame(idx):end)<=half_amp(idx);
            [~, timeDown2halfAmp(idx)] = nanmax(cross_halfAmp);
        end
    end
    metrics.timeUp2halfAmp_nfkb = (timeUp2halfAmp - 1)/12;
    metrics.timeDown2halfAmp_nfkb = (timeDown2halfAmp- 1)/12;

% Half max in early activity, etc
    [max_amp_2h_nfkb, metrics.time2Max_nfkb, metrics.timeUp2halfMax_nfkb, metrics.timeDown2halfMax_nfkb] = early_activity(smoothed(:, StimulationTimePoint:end));          

    function [max_amp, time2Max, timeUp2halfMax, timeDown2halfMax] = early_activity(time_series)
        [max_amp, time2Max] = nanmax(time_series(:,1:25),[],2);
        half_max = max_amp./2;
        cross_halfMax = time_series-half_max>=0;
        [~, timeUp2halfMax] = nanmax(cross_halfMax, [], 2);
        timeDown2halfMax = nan(size(time2Max));
        for index=1:length(time2Max)
            cross_halfMax = time_series(index, time2Max(index):end)<=half_max(index);
            [~, timeDown2halfMax(index)] = nanmax(cross_halfMax);
        end
        timeUp2halfMax = (timeUp2halfMax - 1)/12;
        timeDown2halfMax = (timeDown2halfMax - 1)/12;
        time2Max = (time2Max - 1)/12;
        % if max_amp below time zero value --> timeUp/down will be zero
    end

%}
%% KTR METRICS

if strcmpi(p.Results.IncludeKTR,'on')
%% BASIC KTR METRICS: TIME SERIES, BASELINE, DERIVATIVE, INTEGRAL
% Basic time series. Interpolate over "normal" interval (12 frames per hr) if required
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
    metrics.time_series_ps_ktr = metrics.time_series_ktr(StimulationTimePoint:p.Results.TrimFrame);

% Baseline calculated in filter function into metrics
    metrics.baseline_ktr = info.ktr_baseline;

% Integrated activity
    % includes below 0 activity 
    % start with stimulation time point
    metrics.integrals_ktr = cumtrapz(t(StimulationTimePoint:end),metrics.time_series_ktr(:, StimulationTimePoint:end),2);

% total activity = 'maximum integral'
    metrics.max_integral_ktr = nanmax(metrics.integrals_ktr,[],2);

% integrals above 0 only
    %integrals now are already calculated using only StimulationTimePoint onwards 
    integrals_pos = nan(size(metrics.integrals_ktr));
    integrals_pos(:,1) = metrics.integrals_ktr(:,1);
    for i = 1:size(integrals_pos,1)
        for j = 2:size(metrics.integrals_ktr,2)
            if metrics.integrals_ktr(i,j)>= metrics.integrals_ktr(i,j-1)
             integrals_pos(i,j) =  integrals_pos(i,j-1) + (metrics.integrals_ktr(i,j) - metrics.integrals_ktr(i,j-1));
            elseif metrics.integrals_ktr(i,j)< metrics.integrals_ktr(i,j-1)
             integrals_pos(i,j) =  integrals_pos(i,j-1);
            end
        end
    end

    metrics.integrals_pos_ktr=integrals_pos; 

% total activity above 0= 'maximum integral'
    metrics.max_pos_integral_ktr= nanmax(integrals_pos,[],2);

% time2XMaxActivity (within 1st 8 hours or less)
    endFrame = min(96, size(metrics.integrals_ktr,2));
    P25MaxIntegral = nanmax(integrals_pos(:,1:endFrame),[],2)*0.25; % 25% max integral
    P50MaxIntegral = nanmax(integrals_pos(:,1:endFrame),[],2)*0.5; % half max integral
    P75MaxIntegral = nanmax(integrals_pos(:,1:endFrame),[],2)*0.75; % 75% max integral

    distances25 = abs(integrals_pos(:,1:endFrame)- P25MaxIntegral);
    distances50 = abs(integrals_pos(:,1:endFrame)- P50MaxIntegral);
    distances75 = abs(integrals_pos(:,1:endFrame)- P75MaxIntegral);
    [~, idx25] = nanmin(distances25,[],2);
    [~, idx50] = nanmin(distances50,[],2);
    [~, idx75] = nanmin(distances75,[],2);
    idx25(idx25==1) = NaN;
    idx50(idx50==1) = NaN;
    idx75(idx75==1) = NaN;
    metrics.time2XmaxPosInt_ktr = nan(size(metrics.integrals_pos_ktr, 1), 3);
    metrics.time2XmaxPosInt_ktr(:,1)= (idx25-1)/FramesPerHour;
    metrics.time2XmaxPosInt_ktr(:,2)= (idx50-1)/FramesPerHour;
    metrics.time2XmaxPosInt_ktr(:,3)= (idx75-1)/FramesPerHour;
%%
% Differentiated activity - use central finite difference
    %use StimulationTimePoint as start
    % end after 4h, i.e. 61 TPs (12 baseline, 1 start, 48 2h stim)
    smoothed = smoothdata(fillmissing(metrics.time_series_ktr,'linear', 2,'EndValues','extrap'),2,'sgolay', 6); 
    
    metrics.derivatives_ktr = (smoothed(:,StimulationTimePoint+2:62) - smoothed(:,StimulationTimePoint:60))/(1/6);


%% TRIM EVERYBODY to a common length (of "good" sets, current minimum is roughly 21 hrs)
try
    metrics.time_series_ktr = metrics.time_series_ktr(:,1:p.Results.TrimFrame);
    metrics.time_series_ktr_no_base_ded = metrics.time_series_ktr_no_base_ded(:,1:p.Results.TrimFrame);
    metrics.time_series_ps_ktr = metrics.time_series_ktr(:,StimulationTimePoint:p.Results.TrimFrame);
    metrics.integrals_ktr = metrics.integrals_ktr(:,1:(p.Results.TrimFrame-StimulationTimePoint));
    metrics.integrals_pos_ktr = metrics.integrals_pos_ktr(:,1:(p.Results.TrimFrame-StimulationTimePoint));
%    metrics.derivatives_ktr = metrics.derivatives_ktr(:,1:(p.Results.TrimFrame-2));
    smoothed = smoothed(:,1:p.Results.TrimFrame);
    t = t(1:p.Results.TrimFrame);
catch me
    disp(['Note: vectors too short to cap @ ',num2str(p.Results.TrimFrame),' frames'])
end
    metrics.time_series_smoothed_ktr = smoothed;
%%

% Integrals within half-hour(0-0.5, 0.5-1, 1-1.5), one-hour windows (0-1, 1-2, 2-3) after stimulation and three hour windows (0-3, 1-4, etc) of activity after stimulation
    % automatically only includes time after stimulation (as long as t is adjusted to have 0 at StimulationTimePoint
    max_hr = floor(max(t));
    metrics.intwin1_ktr = nan(size(metrics.time_series_ktr,1),max_hr);
    metrics.intwin3_ktr = nan(size(metrics.time_series_ktr,1),max_hr-2);
    for i = 1:(max_hr)
        win = t>=(i-1) & t<=(i);
        metrics.intwin1_ktr(:,i) = trapz(t(win),metrics.time_series_ktr(:,win),2);
        if i<= (max_hr-2)
            win = t>=(i-1) & t<=(i+2);
            metrics.intwin3_ktr(:,i) = trapz(t(win),metrics.time_series_ktr(:,win),2);
        end
    end
    max_half_hr = floor(max(t)*2)/2;
    metrics.intwin_p5_ktr= nan(size(metrics.time_series_ktr,1),max_half_hr/0.5);
    counter = 1;
    for i = 0.5:0.5:max_half_hr
        win = t>=(i-0.5) & t<=(i); %eg first window goes from 0 to 0.5 h post stim
        metrics.intwin_p5_ktr(:,counter) = trapz(t(win),metrics.time_series_ktr(:,win),2);
        counter = counter + 1;
    end
    %Difference and ratio between hour 0-1 integral and hour 1-2 integral
    metrics.phase_diff1_ktr = metrics.intwin1_ktr(:, 1)-metrics.intwin1_ktr(:, 2);
    metrics.phase_ratio1_ktr = metrics.intwin1_ktr(:, 1)./metrics.intwin1_ktr(:, 2);
    
    %Difference and ratio between hour 0-3 integral and hour 3-6 integral
    metrics.phase_diff3_ktr = metrics.intwin3_ktr(:, 1)-metrics.intwin3_ktr(:, 4);
    metrics.phase_ratio3_ktr = metrics.intwin3_ktr(:, 1)./metrics.intwin3_ktr(:, 4);
%% MAX/MIN metrics
    %include only time after stimulation
    %include only first 4 h
    if size(metrics.time_series_ktr,2)>48+StimulationTimePoint
        metrics.max_amp_ktr = nanmax(metrics.time_series_ktr(:,StimulationTimePoint:48+StimulationTimePoint),[],2);
    end
    
    metrics.max2minDiff_ktr = peak2peak(smoothed(:,StimulationTimePoint:end),2);
%% SIGNAL STATS metrics
    metrics.peak2rms_ktr    =peak2rms(smoothed(:,StimulationTimePoint:end),2);     
    metrics.rms_ktr         =rms(smoothed(:,StimulationTimePoint:end),2);
    metrics.mean_movmad_ktr  =mean( movmad(smoothed(:,StimulationTimePoint:end),3,2),2);
    metrics.mean_movvar_ktr       =mean( movvar(smoothed(:,StimulationTimePoint:end),3,[],2),2);
%% ACTIVITY metrics: compute an off-time for all cells
    %Responder status using sigma threshold and off times
    Wliml = StimulationTimePoint+1; %first/lower time point of window to check for activity
    if size(metrics.time_series_ktr,2)>36+StimulationTimePoint
        Wlimu = StimulationTimePoint + 36; %last/upper time point of window to check for activity, ie check in the first 3 hours after stimulation
    else
        Wlimu = size(metrics.time_series_ktr,2);
    end
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
                metrics.responder_status_ktr(jj,1) = any(block_length>=blockLengthThresh, 2);
                block_length_readout(jj, 1:numel(block_length)) = block_length;
    end
    
    metrics.responders_fraction_ktr = nnz(metrics.responder_status_ktr)/numel(metrics.responder_status_ktr);

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
%% KTR METRICS OF DURATION
% Envelope width: maximum consecutive time above a threshold (envelope must begin within 1st 6 hrs)
    upperThresh = 0.322;
    aux.ktr.thresholds = linspace(0, upperThresh, 5);
    metrics.envelope_ktr = zeros(size(metrics.time_series_ktr,1),length(aux.ktr.thresholds));
    for j = 1:length(aux.ktr.thresholds)
        thresholded = smoothed(:,StimulationTimePoint:end)>aux.ktr.thresholds(j); %consider only TP after stimulation
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

% Duration: Total number of frames above a given threshold
    metrics.duration_ktr = zeros(size(metrics.time_series_ktr,1),length(aux.ktr.thresholds));
    %only include TP from stimualtions onwards
    for i = 1:length(aux.ktr.thresholds)
        metrics.duration_ktr(:,i) = nansum(smoothed(:,StimulationTimePoint:end)>aux.ktr.thresholds(i),2)/FramesPerHour;
    end

%% KTR OSCILLATION METRICS
% Fourier distribution (via FFT) & power
    off_pad = FramesPerHour; % Signal time added to trajectory in  FFT calculation (keeps transients from being recorded as osc.)
    Fs = 1/300;
    depth = max(metrics.off_times_ktr)*FramesPerHour;
    NFFT = 2^nextpow2(depth); % Next power of 2 from chosen depth
    aux.ktr.fft = zeros(size(metrics.time_series_ktr,1),NFFT/2+1);
    aux.ktr.freq = Fs/2*linspace(0,1,NFFT/2+1);
    aux.ktr.power = zeros(size(aux.ktr.fft));
    aux.ktr.power = zeros(size(aux.ktr.fft));
    
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
                   Y = fft(y,NFFT);
                aux.ktr.fft(i,:) = abs(Y(1:NFFT/2+1));
                aux.ktr.power(i,:) = abs(Y(1:NFFT/2+1).^2);
                %SL20200820 Normalization (Plancherel's theorem) to make power peaks comaparable
               aux.ktr.power(i,:) = abs(Y(1:NFFT/2+1).^2)./sum((abs(y)).^2);
            end
        end
    end

% Original Peakfreq metric using aux.power
%{
% Find the point of peak (secondary) power
metrics.peakfreq_ktr = nan(size(aux.ktr.power,1),1);
for i =1:size(metrics.time_series_ktr,1)
    [pks,locs] = globalpeaks(aux.ktr.power(i,:),2);
    % Ensure we're not getting a totally spurious peak
    
    %todo check what to do/how to deal with a number of similarly sized
    %peaks in fft that represent 'oscillating' noise!
    
    %if the two found peaks are not obviously different heights, both peaks kept
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
%}

    % New Peakfreq metric using normalized aux.power
    % Find the point of peak (secondary) power
    metrics.peakfreq_ktr = nan(size(aux.ktr.power,1),1);
    for i =1:size(metrics.time_series_ktr,1)
        [pks,locs] = globalpeaks(aux.ktr.power(i,:),2);
        % Ensure we're not getting a totally spurious peak 
        %if the two found peaks are not obviously different heights, both peaks kept
        if min(pks) < (0.1*max(pks))
            locs(pks==min(pks)) = [];
        elseif min(pks)< 6.5
            locs(pks==min(pks))=[];
        end
        if max(pks)< 6.5
            locs=[];
        end
        %if both peaks kept, use higher frequency peak as peakfreq
        if length(locs)>1
            idx = max(locs(1:2));
            metrics.peakfreq_ktr(i) = 3600*aux.ktr.freq(idx);
        %if locs has 1 value, it takes either frequency corresponding to peak
        %or 3rd frequency (close to 0), whichever is larger, to avoid mixing
        %the low peakfreq determined here becoming indistinguishable from the
        %non-responders with peakfreq = 0
    
        elseif ~isempty(locs)
             metrics.peakfreq_ktr(i) = 3600*aux.ktr.freq(max([locs,3]));
        %if no peaks, peakfreq is 0
        else
            metrics.peakfreq_ktr(i) = 3600*aux.ktr.freq(1);
        end
    end

% Total oscillatory content of particular cells (using thresholds from 0.35 to 0.7 hrs^(-1))

    freq_thresh = aux.ktr.freq( (aux.ktr.freq >= (0.35/3600)) & (aux.ktr.freq <= (0.7/3600)));
%{
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
%}

    %using normalized power
    metrics.oscfrac_ktr = nan(size(aux.ktr.power,1),length(freq_thresh));
    for j = 1:length(freq_thresh)
        for i =1:size(metrics.time_series_ktr,1)
            metrics.oscfrac_ktr(i,j) = nansum(aux.ktr.power(i,((aux.ktr.freq >= freq_thresh(j))&(aux.ktr.power(i,:)>= 6.5)))) /nansum(aux.ktr.power(i,aux.ktr.power(i,:)>= 6.5));
    %        metrics.oscfrac_ktr(i,j) = nansum(aux.ktr.power(i,aux.ktr.freq >= freq_thresh(j))) /nansum(aux.ktr.power(i,:));
            if isnan(metrics.oscfrac_ktr(i,j))
                metrics.oscfrac_ktr(i,j) = 0;
            end
        end
    end

% KTR Oscillation power
    Fs=12; 
    freq_range =[0.35 1];
    n = size(smoothed(:,StimulationTimePoint:end),2); 
    [pwr,fq]=pwelch(smoothed(:,StimulationTimePoint:end)',n,10,256,Fs,'one-sided','psd'); %window size is full length of data? with 10 overlap, 256 is nfft
    fq_ktr =fq';
    %normalize power/psd to 1
    psd_ktr =transpose(pwr./sum(pwr,1));        
    %oscpower, also referred to bandpower
    pwr = transpose(psd_ktr) ; fq = transpose(fq_ktr);% 
    %todo see if freq_range etc need any adjustments esp for KTR
    bp= bandpower(pwr,fq,freq_range, 'psd')';
    metrics.oscpower_ktr=bp;

% 99% oscillation bandwidth 
    metrics.oscbandwidth_ktr =obw(smoothed(:,StimulationTimePoint:end)',Fs)';

%% KTR METRICS OF AMPLITUDE AND TIMING
% 1st peak time/amplitude/prominence/width
    % look for 1st peak in 4h after stimulation
    pk_feats = {'pk1_amp_ktr', 'pk1_time_ktr', 'pk1_width_ktr', 'pk1_prom_ktr'};
    for i=1:length(pk_feats)
        metrics.(pk_feats{i}) = nan(size(metrics.time_series_ktr,1),1);
    end
    
    for i = 1:size(metrics.pk1_time_ktr,1)    
    %   findpeaks(smoothed(i,StimulationTimePoint:(48+StimulationTimePoint)),'NPeaks',8, 'SortStr', 'none', 'WidthReference', 'halfheight', 'MinPeakWidth', 2, 'MinPeakHeight', 0, 'MinPeakProminence',0.5*baseline_stdv_ktr(i),'Annotate', 'extent');
        [pks_ktr, locs_ktr, width_ktr, prom_ktr] = findpeaks(smoothed(i,StimulationTimePoint:(48+StimulationTimePoint)),'NPeaks',8, 'SortStr', 'none', 'WidthReference', 'halfheight', 'MinPeakWidth', 2, 'MinPeakHeight', 0,'MinPeakProminence',0.5*baseline_stdv_ktr(i),'Annotate', 'extent');
       if ~isempty(locs_ktr)
            metrics.pk1_time_ktr(i) = locs_ktr(1);
            metrics.pk1_amp_ktr(i) = pks_ktr(1);
            metrics.pk1_width_ktr(i) = width_ktr(1);
            metrics.pk1_prom_ktr(i) = prom_ktr(1);
       end
    end
    metrics.pk1_time_ktr = ((metrics.pk1_time_ktr)-1)/FramesPerHour;

%% ktr METRICS OF PERC Peak 1 AMP DURATION
% Envelope width: maximum consecutive time above a threshold (envelope must begin within 1st 6 hrs)
    % threshold based on percent max amp for each cell
    upperThresh = 0.8;
    aux.ktr.thresholds = linspace(0.2, upperThresh, 5);
    metrics.envelope_perc_pk1amp_ktr = zeros(size(metrics.time_series_ktr,1),length(aux.ktr.thresholds));
    for j = 1:length(aux.ktr.thresholds)
        thresholded = smoothed(:,StimulationTimePoint:end)>aux.ktr.thresholds(j).*metrics.pk1_amp_ktr;
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
                    if (idx_stop-idx_start) > metrics.envelope_perc_pk1amp_ktr(i,j)
                        metrics.envelope_perc_pk1amp_ktr(i,j) = (idx_stop-idx_start);
                        metrics.starttime_env_perc_pk1amp_ktr(i,j) = idx_start;
                    end
                    curr = idx_stop;
                else
                    break
                end
            end
            if metrics.envelope_perc_pk1amp_ktr(i,j) ==0
                metrics.starttime_env_perc_pk1amp_ktr(i,j) = NaN;
            end
        end  
    end
    metrics.envelope_perc_pk1amp_ktr = metrics.envelope_perc_pk1amp_ktr/FramesPerHour;
    metrics.starttime_env_perc_pk1amp_ktr = metrics.starttime_env_perc_pk1amp_ktr/FramesPerHour;

% Number of frames above a given threshold
    metrics.duration_perc_pk1amp_ktr = zeros(size(metrics.time_series_ktr,1),length(aux.ktr.thresholds));
    for i = 1:length(aux.ktr.thresholds)
        metrics.duration_perc_pk1amp_ktr(:,i) = nansum(smoothed(:,StimulationTimePoint:end)>aux.ktr.thresholds(i).*metrics.pk1_amp_ktr,2)/FramesPerHour;
    end

%% KTR Speed metrics
% Max derivative before pk1 and min derivative within 1 h after pk1, up and down time to half peak amp
    pk1_frame = metrics.pk1_time_ktr *FramesPerHour;
    metrics.max_derivative_pk1_ktr = nan(size(metrics.pk1_time_ktr));
    max_derivative_pk1_frame = nan(size(metrics.pk1_time_ktr));
    metrics.min_derivative_pk1_ktr = nan(size(metrics.pk1_time_ktr));
    min_derivative_pk1_frame = nan(size(metrics.pk1_time_ktr));
    for i = 1:numel(metrics.max_derivative_pk1_ktr )
        if ~isnan(metrics.pk1_time_ktr(i)) &&  pk1_frame(i) <= 36
                [metrics.max_derivative_pk1_ktr(i), max_derivative_pk1_frame(i)] = nanmax(metrics.derivatives_ktr(i,1:pk1_frame(i)),[],2);
                [metrics.min_derivative_pk1_ktr(i), min_derivative_pk1_frame(i)] = nanmin(metrics.derivatives_ktr(i,pk1_frame(i):pk1_frame(i)+12),[],2);
        end
    end
    max_pk1_speed_time_ktr = max_derivative_pk1_frame./FramesPerHour;
    min_pk1_speed_time_ktr = min_derivative_pk1_frame./FramesPerHour;


    metrics.timeUp2halfAmp_ktr = nan(size(metrics.pk1_time_ktr));
    metrics.timeDown2halfAmp_ktr = nan(size(metrics.pk1_time_ktr));
    half_amp = metrics.pk1_amp_ktr./2;
    cross_halfAmp = smoothed(:,StimulationTimePoint:end)-half_amp>=0;
    [~, timeUp2halfAmp] = nanmax(cross_halfAmp, [], 2);
    timeDown2halfAmp = nan(size(pk1_frame));
    
    for idx=1:length(metrics.pk1_time_ktr)
        if ~isnan(metrics.pk1_time_ktr(idx))
            cross_halfAmp = smoothed(idx, StimulationTimePoint+pk1_frame(idx):end)<=half_amp(idx);
            [~, timeDown2halfAmp(idx)] = nanmax(cross_halfAmp);
        end
    end
    metrics.timeUp2halfAmp_ktr = (timeUp2halfAmp- 1)/12;
    metrics.timeDown2halfAmp_ktr = (timeDown2halfAmp - 1)/12;

% Half max in early activity
    [max_amp_2h_ktr, metrics.time2Max_ktr, metrics.timeUp2halfMax_ktr, metrics.timeDown2halfMax_ktr] = early_activity(smoothed(:, StimulationTimePoint:end));          

end
end
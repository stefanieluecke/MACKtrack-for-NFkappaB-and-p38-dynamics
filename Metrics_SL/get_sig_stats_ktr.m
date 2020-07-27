function sig_stats = get_sig_stats_ktr(time_series, StimulationTimePoint,varargin)

%todo adjust to timeframe!!!

%extraFeats =convertCharsToStrings({'power', 'medfreq', 'meanfreq', 'psd', 'noise_est', 'max_pentropy'})';
%is_valid_stat=@(x) all(ismember(x, [string(get_sig_stat_list()); extraFeats])); 
valid_FillMethod  = {'previous', 'next','linear', 'spline', 'pchip'};
isvalidfill =@(x)ismember (x, valid_FillMethod);
p=inputParser;
addRequired(p, 'time_series',@isnumeric);
addRequired(p, 'StimulationTimePoint',@isnumeric);
addParameter(p, 'FqRange',[0.33 1] , @isnumeric);
addParameter(p, 'Fs',12, @isnumeric);
addParameter(p, 'FillMethod','linear', isvalidfill);
%addParameter(p, 'Stats', ["peak2peak", "peak2rms","obw","bandpower"],is_valid_stat);
addParameter(p, 'SmoothMethod',['sgolay'], @istext); %recommended by Ade: 'sgolay', or 'lowess'

parse(p,time_series,StimulationTimePoint, varargin{:});

sig_stats=struct; 

time_series_mod = time_series(:,StimulationTimePoint:end);
time_series_mod=fillmissing(time_series_mod,p.Results.FillMethod, 2, 'EndValues','extrap');
%todo test sgolay and lowess smoothing methods
%todo find out which window size Ade used here
smoothData = smoothdata(time_series_mod,p.Results.SmoothMethod); 

%todo is Fs the same as frames per hour, if so use my parametrization
%todo understand what FqRange does exactly
Fs=p.Results.Fs; freq_range =p.Results.FqRange;

sig_stats.medfreq_ktr      = medfreq(smoothData', Fs,freq_range)';
sig_stats.meanfreq_ktr     = meanfreq(smoothData', Fs, freq_range)'; 
sig_stats.peak2rms_ktr     =peak2rms(smoothData,2);     
sig_stats.rms_ktr          =rms(smoothData')';
sig_stats.peak2peak_ktr    = peak2peak(smoothData,2);
sig_stats.mean_movmad_ktr  =mean( movmad(smoothData,3,2),2);
sig_stats.mean_movstd       =mean( movstd(smoothData,3,[],2),2);
sig_stats.mean_movvar       =mean( movvar(smoothData,3,[],2),2);
  
%psd = power spectral density, also called power in Ade's metrics
    %[~, sig_stats.(stat)]=pwelch(Data',[],[],[],Fs,'one-sided','power');
    %[pwr, fq] = periodogram(smoothData',[], [], Fs, 'power');      
    n = size(time_series_mod,2); 
    [pwr,fq]=pwelch(smoothData',n,10,256,Fs,'one-sided','psd');
%    [pwr,fq]=pwelch(smoothData',n,10,256,Fs,'one-sided',stat);
    sig_stats.fq_ktr    =fq';
    %normalize power/psd to 1
sig_stats.psd_ktr   =transpose(pwr./sum(pwr,1));
            
%oscpower, also referred to bandpower
    pwr = transpose(sig_stats.psd_ktr) ; fq = transpose(sig_stats.fq_ktr);% 
    %todo see if freq_range etc need any adjustments esp for KTR
    bp= bandpower(pwr,fq,freq_range, 'psd')';
sig_stats.oscpower_ktr =bp;

%oscillation frequency
    %find peaks within the frequency range
    bandfilter= @(x) x<= max(freq_range) & x>= min(freq_range);normalize =@(x) x/sum(x);
    ix =bandfilter(sig_stats.fq_ktr);

    %todo Figure out whether peakFun works in teh current context 
    peakFun =@(a) arrayfun(@(j) findpeaks(a(j,:), sig_stats.fq_ktr(ix),...
                    'SortStr', 'descend', 'MinPeakProminence', 0.0055), 1:size(a,1), 'UniformOutput',false);
    [peaks,locs] = peakFun(sig_stats.psd_ktr(:,ix)) ; %peaks = psd, locs = frequency
    freq =zeros(size(peaks)); 
        for j = 1:numel(peaks)
            if numel(peaks{j}) > 1
            % more than one peak within the range, take weighted
            % sum of frequency 
                wgts = normalize(peaks{j}); 
                freq(j) = sum(locs{j}.*wgts); 
            elseif ~isempty(peaks{j})
                freq(j) = locs{j}; 
            end
         end
sig_stats.oscfreq_ktr = freq';
            
            
sig_stats.oscbandwidth_ktr     =obw(smoothData',Fs)';

%max entropy
    max_entropy= zeros(size(smoothData,1),1);    
    time_pts = max_entropy;

        for j =1:numel(max_entropy)
          [sig_entropy, tp]= pentropy(smoothData(j,:), Fs/3600,'FrequencyLimit',freq_range./3600);   
           [max_entropy(j), ix] = max(sig_entropy);
           time_pts(j) = tp(ix)/3600;
        end
sig_stats.max_entropy_ktr = max_entropy; 

sig_stats.noise_est_ktr =  wnoisest(time_series)';

            
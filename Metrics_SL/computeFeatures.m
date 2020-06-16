function [features, metrics]= computeFeatures(id, varargin)

%function that takes experimental ID number and computes features (optionally given 'FeatureList' table) 
%from basic metrics function and additional functions
%(derived from Ade's computeFeatures function, Brooks' nfkbmetrics function, nfkbpeaks function, etc.


%% INPUT PARSING
% Create input parser object, add required params from function input
p = inputParser;
% Required: ID input
valid_id = @(x) assert((isnumeric(x)&&length(x)==1)||exist(x,'file'),...
    'ID input must be spreadsheet ID or full file path');
addRequired(p,'id',valid_id);
% Optional parameters to be passed to metrics function
expectedFlags = {'on','off'};
addParameter(p,'Verbose','on', @(x) any(validatestring(x,expectedFlags)));%checks whether optional name-value argument matches on or off %checks if x matches expectedFlags
valid_conv = @(x) assert(isnumeric(x)&&(x>=0)&&(length(x)==1),...
    'Parameter must be single integer >= 0'); %checks whether parameters below are single integers
addParameter(p,'ConvectionShift',1, valid_conv); %allows adjustment of convection shift (?)
addParameter(p,'MinLifetime',109, @isnumeric); %allows adjustment of minimum lifetime (?)
addParameter(p,'MinSize',90, valid_conv); %allows adjustment of minimum size (?)
addParameter(p,'TrimFrame',157, @isnumeric);
addParameter(p,'StartThreshNFkB',14, valid_conv);  %max allowable starting threshhold (before baseline deduction)to filter out cells with pre-activated NFkB
addParameter (p, 'OnThreshNFkB', 3, @isnumeric); %sigma threshold for determining responders
addParameter(p,'StartThreshKTR',0.9, valid_conv); %max allowable starting threshhold to filter out cells with pre-activated KTR, default is 0.6
addParameter (p, 'OnThreshKTR', 3, @isnumeric); %sigma threshold for determining responders
addParameter(p, 'StimulationTimePoint', 13, @isnumeric)
addParameter(p, 'FramesPerHour', 12, @isnumeric)


% Optional parameters to be passed used in computeFeature function
addParameter(p, 'FeatureListFile', 'C:\Users\stlue\OneDrive\PostDoc UCLA\1 Post Doc UCLA\Matlab analysis\MACKtrack_SL\Metrics_SL\FeatureList.xlsx') %provide file path for Excel table with list of feature to be computed

parse(p,id, varargin{:})

StimulationTimePoint = p.Results.StimulationTimePoint; 
FramesPerHour = p.Results.FramesPerHour;
%%
%get list of features to be calculated
FeatureListFile = p.Results.FeatureListFile;
%FeatureListTable = readtable(FeatureListFile, 'Sheet', 1);
FeatureListTable = readtable(FeatureListFile, 'Sheet', 'FeatureList');
FeatureListTable = readtable(FeatureListFile, 'Sheet', 'ForClassifier');
FeatureList = table2cell(FeatureListTable(:,1));

PeakStatsTable = readtable(FeatureListFile, 'Sheet', 'PeakStats');
PeakStatsList = table2cell(PeakStatsTable(:,1));
SignalStatsTable = readtable(FeatureListFile, 'Sheet', 'SignalStats');
SignalStatsList = table2cell(SignalStatsTable(:,1));


%% call function for basic metrics
[metrics,aux, graph, info, measure] = nfkb_ktr_ratio_metrics(id, 'MinLifetime',p.Results.MinLifetime,...
                            'ConvectionShift',p.Results.ConvectionShift, 'OnThreshNFkB',p.Results.OnThreshNFkB,'OnThreshKTR',p.Results.OnThreshKTR,...
                            'MinSize', p.Results.MinSize,'StartThreshNFkB', p.Results.StartThreshNFkB,'StartThreshKTR', p.Results.StartThreshKTR, 'Verbose', ... 
                            p.Results.Verbose, 'TrimFrame', p.Results.TrimFrame, ...
                            'StimulationTimePoint', p.Results.StimulationTimePoint, 'FramesPerHour', p.Results.FramesPerHour);


%todo call funcitons providing multiple metrics (eg get_peak_stat_list), (if they are requried by Feature List)
% so they're only called once
%todo adjust these inputs to my own fuctions
%get_peak_stat_list
if any(ismember(FeatureList, PeakStatsList))
    peak_stats_nfkb =get_peak_stats_nfkb(metrics.time_series_nfkb, metrics.pk1_amp_nfkb, metrics.pk2_amp_nfkb, 'FramesPerHour', FramesPerHour, 'StimulationTimePoint', StimulationTimePoint);
    peak_stats_ktr =get_peak_stats_ktr(metrics.time_series_ktr, metrics.pk1_amp_ktr, metrics.pk2_amp_ktr, 'FramesPerHour',FramesPerHour, 'StimulationTimePoint', StimulationTimePoint);
end

%sig_stats
if any(ismember(FeatureList, SignalStatsList))
    sig_stats_nfkb =get_sig_stats_nfkb(metrics.time_series_nfkb, StimulationTimePoint);
    sig_stats_ktr =get_sig_stats_ktr(metrics.time_series_ktr, StimulationTimePoint);
end
%get_fold_change
%{
if any(ismember(FeatureList, {'fold_change_nfkb', 'max_fold_change_nfkb'}))
    [fold_change_nfkb, max_fold_change_nfkb] = get_fold_change(metrics.time_series_nfkb_no_base_ded, metrics.baseline_nfkb, StimulationTimePoint);
end
if any(ismember(FeatureList, {'fold_change_ktr', 'max_fold_change_ktr'}))
    [fold_change_ktr, max_fold_change_ktr] = get_fold_change(metrics.time_series_ktr_no_base_ded, metrics.baseline_ktr, StimulationTimePoint);
end
%}

if any(ismember(FeatureList, {'integrals_pos_nfkb','time2HalfMaxPosIntegral_nfkb','max_pos_integral_nfkb',...
                    'integrals_pos_ktr','time2HalfMaxPosIntegral_ktr', 'max_pos_integral_ktr'}))
   % endFrame = min(96+StimulationTimePoint, size(metrics.integrals_nfkb-StimulationTimePoint,2));
     endFrame = min(96+StimulationTimePoint, size(metrics.integrals_nfkb,2)-StimulationTimePoint);
    pos_integral_features_nfkb = get_pos_integrals_nfkb(metrics.integrals_nfkb, FramesPerHour, StimulationTimePoint, endFrame);
    pos_integral_features_ktr = get_pos_integrals_ktr(metrics.integrals_ktr, FramesPerHour, StimulationTimePoint, endFrame);

end



%%
features =struct;

%% 
%get basic metrics, including filtered time_series/trajectories from metrics function
%while some basic metrics, such as time_series, derivatives, integrals consist of all timepoints, the 'feature' should only include measurements
%after stimulation startpoin
for j = 1:length(FeatureList)
    featName = FeatureList{j};
    if ~isfield(features, featName)
        switch featName
%get basic metrics, including filtered time_series/trajectories from metrics function
 %todo integral only above 0 ( in metric function)
            case {'derivatives_nfkb', 'integrals_nfkb', 'time_series_nfkb', 'time_series_nfkb_no_base_ded',...
                    'derivatives_ktr', 'integrals_ktr', 'time_series_ktr', 'time_series_ktr_no_base_ded'}
                features.(featName) = metrics.(featName)(:,StimulationTimePoint:end);
%todo see what to do about intwin metrics (Ade doesn't use them?)
%            case{'intwin1_nfkb','intwin3_nfkb'}
%                features.(featName) = metrics.(featName);
            case{'baseline_nfkb','responder_index_nfkb','off_times_nfkb',...
                    'baseline_ktr','responder_index_ktr','off_times_ktr'}
                features.(featName) = metrics.(featName);
            case{'max_amplitude_nfkb','max_integral_nfkb','max_derivative_nfkb','min_derivative_nfkb',...
                    'max_amplitude_ktr','max_integral_ktr','max_derivative_ktr','min_derivative_ktr'}
                features.(featName) = metrics.(featName);
            case{'peakfreq_nfkb', 'peakfreq_ktr'}
                features.(featName) = metrics.(featName);
            case{'intwin1_nfkb', 'intwin3_nfkb', 'intwin1_ktr','intwin3_ktr'}
                features.(featName) = metrics.(featName);
         %todo osc_frac from basic metrics function, pick threshold or replace entirely with osc_cat?
            case{'oscfrac_nfkb', 'oscfrac_ktr'}
                features.(featName) = metrics.(featName);
            case{'pk1_time_nfkb','pk1_amp_nfkb','pk2_time_nfkb','pk2_amp_nfkb',...
                    'pk1_time_ktr','pk1_amp_ktr','pk2_time_ktr','pk2_amp_ktr'}
                features.(featName) = metrics.(featName);
%todo proper picking of duration and envelope threshold and method,
%different thresholds for KTR And NFkB
            case{'envelope_nfkb','duration_nfkb'}
                features.(featName) = metrics.(featName);
            case{'envelope_sigma_nfkb','duration_sigma_nfkb'}
                features.(featName) = metrics.(featName);
            case{'envelope_ktr','duration_ktr'}
                features.(featName) = metrics.(featName);
            case{'envelope_sigma_ktr','duration_sigma_ktr'}
                features.(featName) = metrics.(featName);
%metrics from Ade's computeFeature functions directly calculated
             case {'pk2_ratio_nfkb'}
                features.(featName)=metrics.pk2_amp_nfkb./metrics.pk1_amp_nfkb;
             case {'pk2_ratio_ktr'}
                features.(featName)=metrics.pk2_amp_ktr./metrics.pk1_amp_ktr;
             case {'median_derivative_nfkb'}
                features.(featName) = nanmedian(metrics.derivatives_nfkb(:,StimulationTimePoint:end),2);
             case {'median_derivative_ktr'}
                features.(featName) = nanmedian(metrics.derivatives_ktr(:,StimulationTimePoint:end),2);
             case {'mean_derivative_nfkb'}
                features.(featName) = nanmean(metrics.derivatives_nfkb(:,StimulationTimePoint:end),2);           
             case {'mean_derivative_ktr'}
                features.(featName) = nanmean(metrics.derivatives_ktr(:,StimulationTimePoint:end),2);
%metrics from Ade's computeFeature functions using additional functions
%todo see if foldchange is useful in my context
            case{'fold_change_nfkb'}
                features.(featName)     = get_fold_change(metrics.time_series_nfkb_no_base_ded, metrics.baseline_nfkb, StimulationTimePoint);
            case{'max_fold_change_nfkb'}
                [~,features.(featName)] = get_fold_change(metrics.time_series_nfkb_no_base_ded, metrics.baseline_nfkb, StimulationTimePoint);
            case{'fold_change_ktr'}
                features.(featName)     = get_fold_change(metrics.time_series_ktr_no_base_ded, metrics.baseline_ktr, StimulationTimePoint);
            case{'max_fold_change_ktr'}
                [~,features.(featName)] = get_fold_change(metrics.time_series_ktr_no_base_ded, metrics.baseline_ktr, StimulationTimePoint);       
%time to half max integral within first 8 hours
            case {'time2HalfMaxIntegral_nfkb'}
                endFrame = min(96+StimulationTimePoint, size(metrics.integrals_nfkb,2));
                features.(featName)     = get_time_to_half_max_integral(metrics.integrals_nfkb(:,1:endFrame), FramesPerHour, StimulationTimePoint);
            case {'time2HalfMaxIntegral_ktr'}
                endFrame = min(96+StimulationTimePoint, size(metrics.integrals_ktr,2));
                features.(featName)     = get_time_to_half_max_integral(metrics.integrals_ktr(:,1:endFrame), FramesPerHour, StimulationTimePoint);
            case{'max_pk1_speed_nfkb'}
                features.(featName)     = get_max_pk1_speed(metrics.pk1_time_nfkb, metrics.derivatives_nfkb, FramesPerHour, StimulationTimePoint);
            case{'max_pk1_speed_ktr'}
                features.(featName)     = get_max_pk1_speed(metrics.pk1_time_ktr, metrics.derivatives_ktr, FramesPerHour, StimulationTimePoint);
            case {'osc_cats_nfkb'}
                 features.(featName)    =  get_osc_cats(metrics.peakfreq_nfkb,metrics.off_times_nfkb,'cutoff_fq', 0.42);
            case {'osc_cats_ktr'}
                 features.(featName)    =  get_osc_cats(metrics.peakfreq_ktr,metrics.off_times_ktr,'cutoff_fq', 0.42);
           
            case{'integrals_pos_nfkb','time2HalfMaxPosIntegral_nfkb','max_pos_integral_nfkb'}
                features.(featName)     = pos_integral_features_nfkb.(featName);
                
            case{'integrals_pos_ktr','time2HalfMaxPosIntegral_ktr', 'max_pos_integral_ktr'}
                features.(featName)     = pos_integral_features_ktr.(featName);
                
                %todo figure out what's up with last_falltime in Ade's code and whether I need it
%           case{'last_falltime'}
%PeakStats            
            case PeakStatsList
                if contains(featName, 'nfkb') 
                    if contains(featName, 'ipt')
                        all_cells= nan(size(peak_stats_nfkb.kept_nfkb, 1),size(peak_stats_nfkb.(featName),2));
                        all_cells(peak_stats_nfkb.kept_nfkb,:)= peak_stats_nfkb.(featName);
                        features.(featName) = all_cells;
                    else
                        features.(featName)=peak_stats_nfkb.(featName);
                    end
                elseif contains(featName, 'ktr')
                    if contains(featName, 'ipt')
                        all_cells= nan(size(peak_stats_ktr.kept_ktr, 1),size(peak_stats_ktr.(featName),2));
                        all_cells(peak_stats_ktr.kept_ktr,:)= peak_stats_ktr.(featName);
                        features.(featName) = all_cells;
                    else
                        features.(featName)=peak_stats_ktr.(featName);
                    end
                end

%SignalStats
            case SignalStatsList
                    if contains(featName, 'nfkb') 
                        features.(featName)=sig_stats_nfkb.(featName);  
                    elseif contains(featName, 'ktr')
                        features.(featName)=sig_stats_ktr.(featName);
                    end
        end
    end
end


%%
%todo Ade distinguishes scalars and vectors --> see whether I need to do that too
%call additional functions to calculate additional metrics

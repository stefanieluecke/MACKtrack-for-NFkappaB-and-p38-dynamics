function [graph, info, measure] = filter_nfkb_ktr_ratio(id,varargin)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% [graph, info, measure] = see_nfkb_ktr(id,graph_flag, verbose_flag)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% SEE_NFKB_KTR is a data processing and visualization script specialized to display
% a nuclear-translocating species (it looks for NFkBdimNuclear and NFkBdimCytoplasm measurements)together with cytoplasmic/nuclear ratio of a second species (e.g. KTR) cell-by-cell.
%
% INPUTS (required):
% id             filename or experiment ID (from Google Spreadsheet specified in "locations.mat")
%
% INPUT PARAMETERS (optional; specify with name-value pairs)
% 'Verbose'         'on' or 'off' - shows non-compacted graphs, e.g.
%                   heatmap including cells to be filtered out(?), default
%                   is off
% 'MinLifetime'     final frame used to filter for long-lived cells, default set to 100
% 'ConvectionShift' Maximum allowable time-shift between different XYs (to correct for poor mixing), default is 1
% 'MinSize'         minimal allowable nuclear area to exclude debris, etc, default set to 90
% 'StartThreshNFkB'     max allowable starting threshhold to filter out cells
%                   with pre-activated NFkB, default is 2
% 'GraphLimitsNFkB' default is [-0.25 8]
% ??'OnThreshNFkB'      default is 1 %? nfkb baseline seems to be defined without this below. Is this used at all..?
% ??'StartThreshKTR'max allowable starting threshhold to filter out cells
%                   with pre-activated KTR, default is 0.6 %currently set
%                   so high no cells are filtered out
%'GraphLimitsKTR'   default is [0 0.35]

% OUTPUTS:  
% graph          primary output structure; must specify
%                   1) filtered/processed data (graph.var) 
%                   2) time vector for all images (graph.t) 
%                   3) XY convection adjustment (graph.shift) 
% info           secondary output structure; must specify
%                   1) Y limits for graphing (info.GraphLimits)
%                   2) parameters from loadID.m (info.parameters) 
% measure         full output structure from loadID
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
%% Create input parser object (checks if functions is provided with correct inputs), add required params from function input
p = inputParser;
% Required: ID input; checks whether ID input is valid (is a numeric array of only
    % one number, or is a structure array (?), or if it is a file)
valid_id = @(x) assert((isnumeric(x)&&length(x)==1)||isstruct(x)||exist(x,'file'),...
    'ID input must be spreadsheet ID or full file path');
addRequired(p,'id',valid_id);

% Optional parameters; specifies allowed optional parameters
expectedFlags = {'on','off'};
addParameter(p,'Verbose','off', @(x) any(validatestring(x,expectedFlags)));%checks whether optional name-value argument matches on or off %checks if x matches expectedFlags
valid_conv = @(x) assert(isnumeric(x)&&(x>=0)&&(length(x)==1),...
    'Parameter must be single integer >= 0'); %checks whether parameters below are single integers
addParameter(p,'ConvectionShift',1, valid_conv); %allows adjustment of convection shift (?)
addParameter(p,'MinLifetime',117, @isnumeric); %allows adjustment of minimum lifetime (?)
addParameter(p,'MinSize',90, valid_conv); %allows adjustment of minimum size (?)
%todo adjust startthresh parameters with new thresholds
addParameter(p,'StartThreshNFkB',10, valid_conv); %max allowable starting threshhold to filter out cells with pre-activated NFkB, default is 2
%addParameter (p, 'OnThreshNFkB', 1, @isnumeric); %? not used in code? --> passed to metric function
%todo adjust OnThresh parameters with new thresholds
addParameter (p, 'GraphLimitsNFkB',[-0.25 8],@isnumeric);
%??? correct start threshholds, etc
addParameter(p,'StartThreshKTR',8, valid_conv); %max allowable starting threshhold to filter out cells with pre-activated KTR, default is 0.6
%addParameter (p, 'OnThreshKTR', 1, @isnumeric);%? not used in code? --> passed to metric function
addParameter (p, 'GraphLimitsKTR',[0 400],@isnumeric);
addParameter(p, 'StimulationTimePoint', 13, @isnumeric)
%2020-05-04 parametrized baseline and background adjustment for NFkB on or off for test function
addParameter(p,'NFkBBaselineDeduction', 'on', @(x) any(validatestring(x,expectedFlags)))
addParameter(p, 'NFkBBackgroundAdjustment', 'on',@(x) any(validatestring(x,expectedFlags)))

% Parse parameters, assign to variables
parse(p,id, varargin{:}) 
if strcmpi(p.Results.Verbose,'on') 
    verbose_flag = 1; %sets verbose_flag to 1 if input parameter Verbose is set to 'on'
else
    verbose_flag = 0;
end

MinLifetime = p.Results.MinLifetime; %pulls number for MinLifetime from input parameters
max_shift = p.Results.ConvectionShift; % Max allowable frame shift in XY-specific correction, takes either default ConvectionShift parameter or the one specified by user

% Set display/filtering parameters
StartThreshNFkB = p.Results.StartThreshNFkB; 
StartThreshKTR = p.Results.StartThreshKTR; 
area_thresh = p.Results.MinSize; % Minimum nuclear area (keeps from including small/junk objects)
StimulationTimePoint = p.Results.StimulationTimePoint;
%% Load AllMeasurements data
[measure, info] = loadID(id);
info.ImageExprNFkB = info.parameters.nfkbdimModule.ImageExpr; %this refers to ImageExpr in nfkbdimModule in parameters of AllMeasurement file
info.ImageExprKTR = info.parameters.ktrModule.ImageExpr;

info.GraphLimitsNFkB = p.Results.GraphLimitsNFkB; % Min/max used in graphing
info.GraphLimitsKTR = p.Results.GraphLimitsKTR; % Min/max used in graphing
info.OnThreshNFkB = p.Results.OnThreshNFkB;
info.OnThreshKTR = p.Results.OnThreshKTR;
baseline_length_nfkb = size(measure.NFkBdim_Nuclear,2); % Endframe for baseline calculation (use entire vector), baseline length is the size of the rows, i.e. number of timepoints 
baseline_length_ktr = size(measure.KTR_ratio1,2); 

%% Filtering
robuststd = @(distr, cutoff) nanstd(distr(distr < (nanmedian(distr)+cutoff*nanstd(distr)))); %standard deviation of ??

% Filtering, part 1 cell fate and cytoplasmic intensity
droprows = []; %creates an empty matrix/array?
droprows = [droprows, sum(isnan(measure.NFkBdim_Nuclear(:,1:4)),2)>2]; % Use only cells existing @ expt start %concatenates a set of 1 or 0 value to droprow matrix (new column?) for each cells depening on whether there are more than 2 NaN values in nuclear NFkB levels within first 4 timepoints
droprows = [droprows, sum(isnan(measure.NFkBdim_Nuclear(:,1:MinLifetime)),2)>3]; % Use only long-lived cells %concatenates a set of 1 or 0 value to droprow matrix (new column?) for each cells depening on whether there are more than 3 NaN values in nuclear NFkB levels within minimum lifetime
%this may need to be adjusted to cyto annu? or keep using full
droprows = [droprows, sum(measure.NFkBdim_Cyto_full(:,1:4)==0,2)>0]; % Very dim cells %concatenates a set of 1 or 0 value to droprow matrix (new column?) for each cells depening on whether there are more than 0 nfkb cytoplasmic values in first four timepoints that are equal to 0
%?? add a filter to remove cells with dim KTR nucleus?
%droprows = [droprows, info.CellData(:,end)]; % Non-edge cells

% NFkB normalization - subtract baseline for each cell (either starting ?% value or 4th percentile of smoothed trajectory)%?whatever is smaller?
nfkb = measure.NFkBdim_Nuclear(:,:); %nfkb is defined as NFkB nuclear measurement from NFkBdim module

%191114 testing differnt versions of baseline

% original baseline calculation by Brooks
%{
nfkb_smooth = nan(size(nfkb)); %NaN array of same size as nfkb is created, for created smoothed trajectory
for i = 1:size(nfkb,1)
    nfkb_smooth(i,~isnan(nfkb(i,:))) = medfilt1(nfkb(i,~isnan(nfkb(i,:))),3); %replaces every element in nfkb_smooth that is not NaN in corresponding nfkb position with a 3rd order median filtered version %?whatever that means...
end
nfkb_min = prctile(nfkb_smooth(:,1:baseline_length_nfkb),5,2); %calculates the 5th percentile along rows of the nfkb smoothed trajectory up to the baseline length, Ade's version uses 2,2 instead of 5,2
nfkb_baseline = nanmin([nanmin(nfkb(:,1:4),[],2),nfkb_min],[],2); %nfkb baseline is defined as minimum of (nfkb_min and the minimum of nfkb at the first four timepoints)(the rows of ?), one per trajectory %?
%?? Parametrize length of baseline?
nfkb = nfkb - repmat(nfkb_baseline,1,size(nfkb,2)); %nfkb is re-defined as nfkb values minus nfkb_baseline array
%}

%
%Testing NFkB baseline, simply deducting mean of first 13 timepoints per cell
if strcmpi(p.Results.NFkBBaselineDeduction,'on')
    %ToDo : Parametrize length of baseline
    nfkb_baseline = nanmean([nfkb(:,1:StimulationTimePoint)],2); %baseline is determined from 1st to 13nth timepoint 
    %? parametrize later on
    nfkb =  nfkb - nfkb_baseline; % ktr activity is defined as baseline - fluorescence measurement
%}
end 

if verbose_flag
    figure, imagesc(nfkb,prctile(nfkb(:),[5,99])),colormap(parula), colorbar %plots raw baseline-subttracted trajectories, using 5th and 99th percentile of nfkb as limits
    title('All (baseline-subtracted) NFkB trajectories')
end

if strcmpi(p.Results.NFkBBackgroundAdjustment,'on')
    nfkb = nfkb/mean(info.parameters.adj_distr_NFkBdim(2,:)); %nfkb is re-defined of nfkb divided by mean of second row of adj_distr %? but I don't know what that does...
    %? why is this done after baseline subtractions, etc?
end   
    
    
ktr = measure.KTR_ratio1(:,:); %ktr is defined as KTR cytoplasmic/nuclear ratio from ktr module

%KTR baseline deduction --> baseline deduction is appropriate only for some
%metrics --> pass along both
ktr_baseline = nanmean([ktr(:,1:StimulationTimePoint)],2); %baseline is determined from 1st to 13nth timepoint 
    %? parametrize later on
ktr_baseline_deducted =  ktr - ktr_baseline; % ktr activity is defined as baseline - fluorescence measurement
%}

%test older baseline 
%{
ktr_smooth = nan(size(ktr)); %NaN array of same size as ktr is created, for created smoothed trajectory
for i = 1:size(ktr,1)
    ktr_smooth(i,~isnan(ktr(i,:))) = medfilt1(ktr(i,~isnan(ktr(i,:))),3); %replaces every element in ktr_smooth that is not NaN in corresponding ktr position with a 3rd order median filtered version %?whatever that means...
end
ktr_min = prctile(ktr_smooth(:,1:baseline_length_ktr),5,2); %calculates the 5th percentile along rows of the ktr smoothed trajectory up to the baseline length, Ade's nfkb version uses 2,2 instead of 5,2

ktr_baseline = nanmin([nanmin(ktr(:,1:4),[],2),ktr_min],[],2); %ktr baseline is defined as minimum of (ktr_min and the minimum of ktr at the first four timepoints)(the rows of ?) %? what is this for?
ktr = ktr - repmat(ktr_baseline,1,size(ktr,2)); %ktr is re-defined as ktr values minus ktr_baseline array
%}

if verbose_flag
    figure, imagesc(ktr,prctile(ktr(:),[5,99])),colormap(parula), colorbar %plots raw baseline-subttracted trajectories, using 5th and 99th percentile of ktr ratio as limits
    title('All (baseline-subtracted) KTR trajectories')
end


%? testing ktr annulus module ratio, leave this out at first
%Todo: This division is likely not necessary appropriate for KTR, because
%ratio is being used --> but test variation in unstimulated points over
%different exp to be sure

%ktr = ktr/mean(info.parameters.adj_distr_ktr(2,:)); %ktr is re-defined of ktr divided by mean of second row of adj_distr_ktr from modebalance in loadID to adjust numbers (fluorscence intensitites) to match between different days
    %? double check this with Ade
    %? Why is this placed here and not before baseline corrections?

% Filtering, part 2: eliminate outlier cells (based on mean value)
nfkb_lvl = reshape(nfkb(max(droprows,[],2) == 0,:),[1 numel(nfkb(max(droprows,[],2) == 0,:))]); %sets nfkb level as nfkb measurements of cells (full trajectories) to be kept reshaped to a 1xtotal number of array elements matrix 
droprows =  [droprows, (nanmean(abs(nfkb-nanmean(nfkb_lvl)),2)./nanstd(nfkb_lvl))>=3]; %removes cells with mean levels larger than 3x standard deviation (mean of all values subtracted from each nfkb element, absolute value of that, mean across each trajectory, each element dividided by standard deviation of nfkb level, check if larger or equal than 3, add to droprows as a 1/0 column
droprows =  [droprows, (nanmean(abs(nfkb-nanmean(nfkb_lvl)),2)./nanstd(nfkb_lvl))>=1.7]; %removes cells with mean levels larger than 1.7x standard deviation 
% Filtering, part 2: eliminate outlier cells (based on mean value)
ktr_lvl = reshape(ktr(max(droprows,[],2) == 0,:),[1 numel(ktr(max(droprows,[],2) == 0,:))]); %sets ktr level as ktr ratio measurements of cells (full trajectories) to be kept reshaped to a 1xtotal number of array elements matrix 
%? test whether these threshold are appropriate for ktr and where they come from
droprows =  [droprows, (nanmean(abs(ktr-nanmean(ktr_lvl)),2)./nanstd(ktr_lvl))>=3]; %removes cells with mean levels larger than 3x standard deviation (mean of all values subtracted from each ktr element, absolute value of that, mean across each trajectory, each element dividided by standard deviation of ktr level, check if larger or equal than 3, add to droprows as a 1/0 column
droprows =  [droprows, (nanmean(abs(ktr-nanmean(ktr_lvl)),2)./nanstd(ktr_lvl))>=1.7]; %removes cells with mean levels larger than 1.7x standard deviation 

% Filtering, part 3: nuclear stain intensity and starting NFkB value
keep = max(droprows,[],2) == 0; %index of cells to keep are those rows of the droprows vector where no columns show a 1
start_lvl_nfkb = nanmin(nfkb(keep,1:3),[],2); % this seems to be used to show nfkb start level, but is not used to filter, calculates start level for each cell (those kept) based on the minimum of the first three timepoints (why the minimum and not the average?)
    %? add in parametrized starting time point?
start_lvl_ktr = nanmin(ktr(keep,1:3),[],2); % this seems to be used to show ktr start level, but is not used to filter, calculates start level for each cell (those kept) based on the minimum of the first three timepoints (why the minimum and not the average?)
%? add in parametrized starting time point?

%
nuc_lvl = nanmedian(measure.MeanNuc1(keep,1:31),2); %defines nuc_lvl as median nuclear intensity (this is DNA staining?) in first 31 timepoints, for each column/timepoint(?) 
nuc_thresh = nanmedian(nuc_lvl)+2.5*robuststd(nuc_lvl(:),2); %defines the acceptable threshold for nuclear intensity as median of nuclear levels + 2.5x fancy std defined above ???

droprows =  [droprows, prctile(nfkb(:,1:8),18.75,2) > StartThreshNFkB]; %removes cells for which the 19th percentile of the first 8 timepoints is larger than start threshold %/why?
%droprows =  [droprows, prctile(ktr(:,1:8),18.75,2) < StartThreshKTR]; %use this when using non-inverted KTR %removes cells for which the 19th percentile of the first 8 timepoints is larger than start threshold %/why?
droprows =  [droprows, prctile(ktr(:,1:8),18.75,2) > StartThreshKTR]; %removes cells for which the 19th percentile of the first 8 timepoints is larger than start threshold %/why?
droprows =  [droprows, nanmedian(measure.MeanNuc1(:,1:31),2) > nuc_thresh];%removes cells with median nuclear intensity above nuclear threshold in first 31 timepoints
droprows =  [droprows, nanmedian(measure.Area,2) < area_thresh]; %removes cells with median areasmaller than the area threshold (to remove small particles)
info.dropped = droprows; %added from Ade's version 06-24-19

% Show some filter information
if verbose_flag
    filter_str = {'didn''t exist @ start', 'short-lived cells', 'very dim NFkB',...
        'NFkB outliers [mean val >3*std]','extreme NFkB val [mean>1.7*std]','extreme KTR val [mean>3*std]', 'KTR outliers [mean val >1.7*std]', 'NFkB active @ start','KTR active @ start', 'high nuclear stain','low area'};
    disp(['INITIAL: ', num2str(size(droprows,1)),' cells'])
    for i = 1:size(droprows,2)
        if i ==1
            num_dropped = sum(droprows(:,i)==1);
        else
            num_dropped = sum( (max(droprows(:,1:i-1),[],2)==0) & (droprows(:,i)==1));
        end
        disp(['Filter #', num2str(i), ' (',filter_str{i},') - ',num2str(num_dropped), ' cells dropped']) 
    end
    disp(['FINAL: ', num2str(sum(max(droprows,[],2) == 0)),' cells'])
%? later, add verbose small multiple graphs if desired
end

info.keep = max(droprows,[],2) == 0;
nfkb = nfkb(info.keep,:); %nfkb is redefined as only the not-filtered-out cells (?)
ktr = ktr(info.keep,:); %ktr is redefined as only the not-filtered-out cells (?)
ktr_baseline_deducted = ktr_baseline_deducted(info.keep,:);
%% Initialize outputs, do final corrections
graph.celldata = info.CellData(info.keep,:); %graph.celldata is set to include celldata from all the non-filtered out cells

%? Add this in later, think about whether it makes sense to have this separately
% Correct for XY positions that activate late
%[graph.var_nfkb, shift_xy_nfkb] = alignTrajectories(nfkb, graph.celldata, 60, max_shift); %calls subfunction alignTrajectories for nfkb, using celldata, 60 as a window to calculate pairwise distances (???), and maximum shift user determined (or default 1 in this function; subfunction sets default to 3)
%[graph.var_ktr, shift_xy_ktr] = alignTrajectories(ktr, graph.celldata, 60, max_shift); %calls subfunction alignTrajectories for ktr, using celldata, 60 as a window to calculate pairwise distances (???), and maximum shift user determined (or default 1 in this function; subfunction sets default to 3)
%temporary definition of graph.var:
graph.var_nfkb = nfkb;
graph.var_ktr = ktr;
graph.var_ktr_baseline_deducted = ktr_baseline_deducted;

%displays the xy-shift done by alignTrajectories for each xy position
%removed because alignTrajectories removed temporarily
%if verbose_flag
%   for i = 1:length(shift_xy)
%        disp(['xy ',num2str(i),' shift : ',num2str(shift_xy(i))])
%   end
%end

%I do not understand this part yet?
%todo removed temporarily, add back in
%graph.shift = shift_xy; %gets the xy shift from alignTrajectories function
%graph.t = 0:(1/info.parameters.FramesPerHour):48;
graph.t = ((-StimulationTimePoint+1)/info.parameters.FramesPerHour):(1/info.parameters.FramesPerHour):48; %creates a time axis vector for the graph from 0 to 48 in steps of 1/FramesperHour (12) %?why
graph.t = graph.t(1:min([length(graph.t),size(graph.var_nfkb,2)]));%time axis vector shortened to number of timepoints in data (if shorter)
%graph.opt_nfkb = maketicks(graph.t,info.GraphLimitsNFkB,0); %calls function to add tick labels to time frame (?only partially understand what that does)
%graph.opt_ktr = maketicks(graph.t,info.GraphLimitsKTR,0); %calls function to add tick labels to time frame (?only partially understand what that does)
%graph.opt_nfkb.Name = 'NFkB activation'; 
%graph.opt_ktr.Name = 'Kinase activation'; 

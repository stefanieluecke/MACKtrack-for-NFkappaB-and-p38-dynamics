function [graph, info, measure] = see_nfkb_native(id,varargin)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% [graph, info, measure] = see_nfkb_native(id,graph_flag, verbose_flag)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% SEE_NFKB_NATIVE is a data processing.and visualization script specialized to handle
% a nuclear-translocating species (it looks for NFkBdimNuclear and NFkBdimCytoplasm measurements).
%
% INPUTS (required):
% id             filename or experiment ID (from Google Spreadsheet specified in "locations.mat")
%
% INPUT PARAMETERS (optional; specify with name-value pairs)
% 'Display'         'on' or 'off' - show graphs (default: process data
%                    only; no graphs), default is off
% 'Verbose'         'on' or 'off' - shows non-compacted graphs, e.g.
%                   heatmap including cells to be filtered out(?), default
%                   is off
% 'MinLifetime'     final frame used to filter for long-lived cells, default set to 100
% 'ConvectionShift' Maximum allowable time-shift between different XYs (to correct for poor mixing), default is 1
% 'MinSize'         minimal allowable nuclear area to exclude debris, etc, default set to 90
% 'StartThresh'     max allowable starting threshhold to filter out cells
%                   with pre-activated NFkB, default is 2
% 'GraphLimits'     default is [-0.25 10]
% 'Baseline'        default is 1 %? nfkb baseline seems to be defined without this below. Is this used at all..?
% 'Source'          default is 'ade' %used in Ade version as input for loadID function, thus unneccessary, unless I add more input arguments when calling loadID

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
addParameter(p,'Display','off', @(x) any(validatestring(x,expectedFlags)));%checks whether optional name-value argument matches on or off %checks if x matches expectedFlags
addParameter(p,'Verbose','off', @(x) any(validatestring(x,expectedFlags)));%checks whether optional name-value argument matches on or off %checks if x matches expectedFlags
valid_conv = @(x) assert(isnumeric(x)&&(x>=0)&&(length(x)==1),...
    'Parameter must be single integer >= 0'); %checks whether parameters below are single integers
addParameter(p,'ConvectionShift',1, valid_conv); %allows adjustment of convection shift (?)
addParameter(p,'MinLifetime',100, @isnumeric); %allows adjustment of minimum lifetime (?)
addParameter(p,'MinSize',90, valid_conv); %allows adjustment of minimum size (?)
%06-20-19 addition of input parameter to easily change starting activity threshold
addParameter(p,'StartThresh',2, valid_conv); %max allowable starting threshhold to filter out cells with pre-activated NFkB, default is 2
%06-21-19 added optional parameters from Ade
addParameter (p, 'Baseline', 1, @isnumeric);
addParameter (p, 'GraphLimits',[-0.25 10],@isnumeric);
addParameter (p, 'Source', 'ade'); %? see whether I need to include this


% Parse parameters, assign to variables
parse(p,id, varargin{:}) %? I don't understand what 'parse' really does
if strcmpi(p.Results.Verbose,'on') 
    verbose_flag = 1; %sets verbose_flag to 1 if input parameter Verbose is set to 'on'
else
    verbose_flag = 0;
end
if strcmpi(p.Results.Display,'on') %sets graph_flag to 1 if input parameter Graph is set to 'on'
    graph_flag = 1; 
else
    graph_flag = 0;
end
MinLifetime = p.Results.MinLifetime; %pulls number for MinLifetime from input parameters
max_shift = p.Results.ConvectionShift; % Max allowable frame shift in XY-specific correction, takes either default ConvectionShift parameter or the one specified by user

% Set display/filtering parameters
% 06-19-19 test for less filtering 
%start_thresh = 2; % Maximal allowable start level above baseline
%start_thresh = 50; % Maximal allowable start level above baseline
start_thresh = p.Results.StartThresh; % 06-20-19 add StartThresh as input parameter
area_thresh = p.Results.MinSize; % Minimum nuclear area (keeps from including small/junk objects)

%% Load AllMeasurements data
[measure, info] = loadID(id);
info.ImageExpr = info.parameters.nfkbdimModule.ImageExpr; %this refers to ImageExpr in nfkbdimModule in parameters of AllMeasurement file

%06-19-19 test for show full activation by changing upper graph limit
%info.GraphLimits = [-0.25 8]; % Min/max used in graphing, replaced 06-21-19
info.GraphLimits = p.Results.GraphLimits; % Min/max used in graphing
info.Baseline = p.Results.Baseline;
%test 06-26-19
% dendro = 1;
dendro = 0;
colors = setcolors;
baseline_length = size(measure.NFkBdim_Nuclear,2); % Endframe for baseline calculation (use entire vector), baseline length is the size of the rows, i.e. number of timepoints 

%{
% Look for locations .mat in default MACKtrack directory, or in current directory
    % Ade's version doesn't include this whole load section --> unnecessary?
home_folder = mfilename('fullpath');
slash_idx = strfind(home_folder,filesep);
if exist([home_folder(1:slash_idx(end-1)),'locations.mat'],'file')
    load([home_folder(1:slash_idx(end-1)), 'locations.mat'],'-mat')
elseif exist([home_folder(1:slash_idx(end)),'locations.mat'],'file')
    load([home_folder(1:slash_idx(end-1)), 'locations.mat'],'-mat')
end
%}
%{
% Experiment-specific visualization settings/tweaks (set by spreadsheet URL)
% BT's experiments
if isnumeric(id)
    if  ~isempty(strfind(locations.spreadsheet,'2PACX-1vSmFCMOlGYOUzq'))
        % a) Heterozygous cell experiments
        if (id <= 270) || ismember(id,[370:379, 384:391, 395, 396])
            start_thresh = 1.5;
            info.GraphLimits = [-0.25 5.5];
        end

        % b) bad light guide (4 experiments from same day)
        if ismember(id,315:318)
            info.GraphLimits = [-0.2 4];
        end
        % d) 100uM CpG - delayed stimulation
        if id==283
            measure.NFkBdimNuclear = measure.NFkBdimNuclear(:,4:end);
            measure.NFkBdimCytoplasm = measure.NFkBdimNuclear(:,4:end);
            disp('Adjusted start point for this CpG expmt')
        end

        % e) prestimulated sets; don't filter pre-activated cells
        if ismember(id,[267:270, 323:324, 337:339,342:343, 356:360])
            start_thresh = 10;
        end

        % f) 09/23/2015 Missed an early timepoint; don't allow back-shifting
        if ismember(id,337:343) 
            max_shift = 0;
        end
        % g) IkBa KO set: poor nuclear staining (weak nuclei predominate): use eroded version of nuclear NFkB
        if ismember(id, 366:379)
            measure.NFkBdimNuclear = measure.NFkBdimNuclear_erode;
            disp('(Using eroded nuclei)')
        end
    % AA's experiments    
    elseif ~isempty(strfind(locations.spreadsheet,'1s51cOI7NvLIOEpEhZWjJKsPkCOE5Qz39m4tHa9nJ7ok'))
        % a) early experiments; heterozygous cells
        if id <= 60
            start_thresh = 1.5;
            info.GraphLimits = [-0.25 6];
        else
            start_thresh = 1.8;
            info.GraphLimits = [-0.25 5.5];
            info.baseline = 0.75;
        end
    end
end
%}
%% Filtering
robuststd = @(distr, cutoff) nanstd(distr(distr < (nanmedian(distr)+cutoff*nanstd(distr)))); %standard deviation of ??

% Filtering, part 1 cell fate and cytoplasmic intensity
droprows = []; %creates an empty matrix/array?
droprows = [droprows, sum(isnan(measure.NFkBdim_Nuclear(:,1:4)),2)>2]; % Use only cells existing @ expt start %concatenates a set of 1 or 0 value to droprow matrix (new column?) for each cells depening on whether there are more than 2 NaN values in nuclear NFkB levels within first 4 timepoints
droprows = [droprows, sum(isnan(measure.NFkBdim_Nuclear(:,1:MinLifetime)),2)>3]; % Use only long-lived cells %concatenates a set of 1 or 0 value to droprow matrix (new column?) for each cells depening on whether there are more than 3 NaN values in nuclear NFkB levels within minimum lifetime
droprows = [droprows, sum(measure.NFkBdim_Cyto_full(:,1:4)==0,2)>0]; % Very dim cells %concatenates a set of 1 or 0 value to droprow matrix (new column?) for each cells depening on whether there are more than 0 nfkb cytoplasmic values in first four timepoints that are equal to 0
%droprows = [droprows, info.CellData(:,end)]; % Non-edge cells

% NFkB normalization - subtract baseline for each cell (either starting ?% value or 4th percentile of smoothed trajectory)%?whatever is smaller?
nfkb = measure.NFkBdim_Nuclear(:,:); %nfkb is defined as NFkB nuclear measurement from NFkBdim module
nfkb_smooth = nan(size(nfkb)); %NaN array of same size as nfkb is created, for created smoothed trajectory
for i = 1:size(nfkb,1)
    nfkb_smooth(i,~isnan(nfkb(i,:))) = medfilt1(nfkb(i,~isnan(nfkb(i,:))),3); %replaces every element in nfkb_smooth that is not NaN in corresponding nfkb position with a 3rd order median filtered version %?whatever that means...
end
nfkb_min = prctile(nfkb_smooth(:,1:baseline_length),5,2); %calculates the 5th percentile along rows of the nfkb smoothed trajectory up to the baseline length, Ade's version uses 2,2 instead of 5,2

nfkb_baseline = nanmin([nanmin(nfkb(:,1:4),[],2),nfkb_min],[],2); %nfkb baseline is defined as minimum of (nfkb_min and the minimum of nfkb at the first four timepoints)(the rows of ?), one per trajectory %?
nfkb = nfkb - repmat(nfkb_baseline,1,size(nfkb,2)); %nfkb is re-defined as nfkb values minus nfkb_baseline array
if verbose_flag
    figure, imagesc(nfkb,prctile(nfkb(:),[5,99])),colormap(parula), colorbar %plots raw baseline-subttracted trajectories, using 5th and 99th percentile of nfkb as limits
    title('All (baseline-subtracted) trajectories')
end
nfkb = nfkb/mean(info.parameters.adj_distr_NFkBdim(2,:)); %nfkb is re-defined of nfkb divided by mean of second row of adj_distr %? but I don't know what that does...
%allow image to image comparison

% Filtering, part 2: eliminate outlier cells (based on mean value)
nfkb_lvl = reshape(nfkb(max(droprows,[],2) == 0,:),[1 numel(nfkb(max(droprows,[],2) == 0,:))]); %sets nfkb level as nfkb measurements of cells (full trajectories) to be kept reshaped to a 1xtotal number of array elements matrix 
droprows =  [droprows, (nanmean(abs(nfkb-nanmean(nfkb_lvl)),2)./nanstd(nfkb_lvl))>=3]; %removes cells with mean levels larger than 3x standard deviation (mean of all values subtracted from each nfkb element, absolute value of that, mean across each trajectory, each element dividided by standard deviation of nfkb level, check if larger or equal than 3, add to droprows as a 1/0 column
droprows =  [droprows, (nanmean(abs(nfkb-nanmean(nfkb_lvl)),2)./nanstd(nfkb_lvl))>=1.7]; %removes cells with mean levels larger than 1.7x standard deviation 

% Filtering, part 3: nuclear stain intensity and starting NFkB value
keep = max(droprows,[],2) == 0; %index of cells to keep are those rows of the droprows vector where no columns show a 1
start_lvl = nanmin(nfkb(keep,1:3),[],2); % this seems to be used to show nfkb start level, but is not used to filter, calculates start level for each cell (those kept) based on the minimum of the first three timepoints (why the minimum and not the average?)

% test :4/25/2019
measure.MeanIntensityNuc = measure.MeanNuc1;

%
nuc_lvl = nanmedian(measure.MeanIntensityNuc(keep,1:31),2); %defines nuc_lvl as median nuclear intensity (this is DNA staining?) in first 31 timepoints, for each column/timepoint(?) 
nuc_thresh = nanmedian(nuc_lvl)+2.5*robuststd(nuc_lvl(:),2); %defines the acceptable threshold for nuclear intensity as median of nuclear levels + 2.5x fancy std defined above ???

droprows =  [droprows, prctile(nfkb(:,1:8),18.75,2) > start_thresh]; %removes cells for which the 19th percentile of the first 8 timepoints is larger than start threshold %/why?
droprows =  [droprows, nanmedian(measure.MeanIntensityNuc(:,1:31),2) > nuc_thresh];%removes cells with median nuclear intensity above nuclear threshold in first 31 timepoints
droprows =  [droprows, nanmedian(measure.Area,2) < area_thresh]; %removes cells with median areasmaller than the area threshold (to remove small particles)
info.dropped = droprows; %added from Ade's version 06-24-19

% Show some filter information
if verbose_flag
    filter_str = {'didn''t exist @ start', 'short-lived cells', 'NFkB<background',...
        'outliers [mean val >3*std]','extreme val [mean>1.7*std]', 'active @ start', 'high nuclear stain','low area'};
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

    %This part creates the small multiples line graphs of individual
    %trajectories ordered by starting value, by nuclear stain level, and by median area
    %current (06/25/19)ranksmult function plots 15x8 subplots with
    %ascending ordering (evenly spaced over entire distribution (ie representing everything from lowest to highest)
    % The title problem is not fully solved, the title is in the wrong place, fix later
    ranksmult(nfkb(keep,:),start_lvl); %ranksmult is function stored under utilities, uses nfkb matrix and the start level for each cell as input 
    %h = suptitle(['x = NFkB activation @ start. Threshold = ',num2str(start_thresh)]);%can't find 'subtitle' in Matlab documentation, but occurs both in Ade's and Brook's version...
    h = sgtitle(['x = NFkB activation at start. Threshold = ',num2str(start_thresh)]);
    set(h,'FontSize',14)
    ranksmult(nfkb(keep,:),nuc_lvl);
    %h = suptitle(['x = Nuclear stain level. Threshold = ',num2str(nuc_thresh)]);
    h = sgtitle(['x = Nuclear stain level. Threshold = ',num2str(nuc_thresh)]);
    set(h,'FontSize',14)
    ranksmult(nfkb(keep,:),nanmedian(measure.Area(keep,:),2));
    %h = suptitle(['x = Median area. Threshold = ',num2str(area_thresh)]);
    h = sgtitle(['x = Median area. Threshold = ',num2str(area_thresh)]);
    set(h,'FontSize',14)

end

info.keep = max(droprows,[],2) == 0;
nfkb = nfkb(info.keep,:); %nfkb is redefined as only the not-filtered-out cells (?)

%% Initialize outputs, do final corrections
graph.celldata = info.CellData(info.keep,:); %graph.celldata is set to include celldata from all the non-filtered out cells

% Correct for XY positions that activate late
[graph.var, shift_xy] = alignTrajectories(nfkb, graph.celldata, 60, max_shift); %calls subfunction alignTrajectories for nfkb, using celldata, 60 as a window to calculate pairwise distances (???), and maximum shift user determined (or default 1 in this function; subfunction sets default to 3)
% this part doesn't run ever, since dendro is set to 0 above and no comparison occurs after 'if'???
% in Ade's version the if dendro part is commented out, but what is the ordering happening after else?
% when the else part is commented out, there's no graph.order, so the code doesn't run --> how does Ade's version run?
if dendro %if dendro is 1, not 0 (maybe make this input parameter?)
    [graph.order, graph.dendro.links] = hierarchial(graph.var(:,1:min([size(graph.var,2),150])),0); %calls function to make order for dendrogram ?don't know what this is for
else
    %sort by total activity over entire time??
    [~,graph.order] = sort(nansum(graph.var(:,1:min([size(graph.var,2),150])),2),'descend'); %sorts the sum of all rows from shift-adjusted trajectories (from timepoint 1 to at most 150 (or amount of timepoints in shift-adjested trajectories) in descending order and creates an index called graph.order from this
end
%displays the xy-shift done by alignTrajectories for each xy position
if verbose_flag
    for i = 1:length(shift_xy)
        disp(['xy ',num2str(i),' shift : ',num2str(shift_xy(i))])
    end
end

%I do not understand this part yet?
graph.shift = shift_xy; %gets the xy shift from alignTrajectories function
graph.t = 0:(1/info.parameters.FramesPerHour):48; %creates a time axis vector for the graph from 0 to 48 in steps of 1/FramesperHour (12) %?why
graph.t = graph.t(1:min([length(graph.t),size(graph.var,2)]));%time axis vector shortened to number of timepoints in data (if shorter)
graph.opt = maketicks(graph.t,info.GraphLimits,0); %calls function to add tick labels to time frame (?only partially understand what that does)
graph.opt.Name = 'NFkB Activation'; 


%% Graphing
if graph_flag
    % Colormap stack: creates heatmap of filtered trajectories using order determined above (graph.order)
    figs.a = figure('name','ColormapStack');
    set(figs.a,'Position', [500 700 1200 600])
    colormapStack(graph.var(graph.order,:),graph.celldata(graph.order,:), graph.opt); %calls subfunction that makes figure    
    
    % Hierarchial linkage
    if dendro
        graph.dendro.label = cell(size(graph.var,1),1);
        graph.dendro.label(:) = {''};
        figs.b = figure('name','EuclideanDist');
        graph.dendro.lines = dendrogram(graph.dendro.links,0,'orientation','left','labels',graph.dendro.label);
        set(figs.b,'Position', [555   743   300   535])
        set(get(graph.dendro.lines(1),'Parent'),'LooseInset',get(get(graph.dendro.lines(1),'Parent'),'TightInset')+ [0 0 0 0])
        set(graph.dendro.lines,'Color',[0.3 0.3 0.3])
    end

    % Line plot (mean+/-std)
    graph.line.top = nanmean(graph.var) + nanstd(graph.var); %creates upper bound
    graph.line.bot = nanmean(graph.var) - nanstd(graph.var); %creates lower bound
    figs.c = figure('name','MeanTrajectory'); %creates figure and names it
    fill([graph.t,graph.t(end:-1:1)],[graph.line.top,graph.line.bot(end:-1:1)],colors.light_blue) %creates filled polygon determined by time axis and upper and lower bound
    hold on
    graph.line.main = plot(graph.t, nanmean(graph.var)); %adds the mean to the figure
    set(graph.line.main,'Color',[.25 .25 .25],'LineWidth',2); %determined properties of the mean line
    hold off
    set(gca,'XTick',graph.opt.TimeTicks,'YTick',graph.opt.MeasurementTicks) %sets axis tick properites
    set(gca,'YTickLabel',graph.opt.MeasurementTickLabels,'TickLength',[0.005 0.005]) %sets tick label properties
    set(figs.c,'Position', [555   743   1145   300]) %positions figure within window
    axis([min(graph.opt.Times) max(graph.opt.Times) graph.opt.MeasurementBounds]) %sets axis to determined timescale and determined measurement bounds, derived from maketicks function

    % Small multiples line graph
    figs.d  = figure('name','smallmultiples');
    set(figs.d,'Position',[500, 350, 876, 1000]);
    graph.smult.order = randperm(size(graph.var,1),min([60 size(graph.var,1)])); %order the plots randomly (vector of max 60 or number of cells random integers between 1 and the number of cells)
    graph.smult.h = tight_subplot(15,4); %creates 15x4 small subplots using function
    xpos = max(graph.opt.Times)-0.48*(max(graph.opt.Times)-min(graph.opt.Times)); %determines the positioning of text within each small graph
    ypos =  max(graph.opt.MeasurementBounds) - 0.26*diff(graph.opt.MeasurementBounds); %determines the positioning of text within each small graph
    for i =1:length(graph.smult.order)
        %plot(graph.smult.h(i),graph.t,graph.var(graph.smult.order(i),:),'Color',colors.grays{3}, 'LineWidth',2)
        plot(graph.smult.h(i),graph.t,graph.var(graph.smult.order(i),:),'Color',colors.red, 'LineWidth',2)
        set(graph.smult.h(i),'XLim',[min(graph.opt.Times)-(range(graph.opt.Times)*0.02) max(graph.opt.Times)],...
            'YLim',graph.opt.MeasurementBounds,'XTickLabel',{},'YTickLabel',{})
        text(xpos,ypos,['XY ',num2str(graph.celldata(graph.smult.order(i),1)),...
            ', cell ',num2str(graph.celldata(graph.smult.order(i),2))],'Parent',graph.smult.h(i))
            
    end
end
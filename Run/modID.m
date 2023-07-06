function [] = modID(varargin)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% [] = runID(varargin)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% RUNID queues and tracks multiple experimental sets from the "Scope Runs" Google Doc -  
% choose sets by ID number.
%
% INPUT:
% varargin     ID# of sets to track (1st column of spreadsheet) - single values or vector
%
% Example Usage: runID(1,4:16,27) will sequentially track sets 1, 4 to 16, and 27.
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

% Specify save directory and load spreadsheet URL
home_folder = mfilename('fullpath');
slash_idx = strfind(home_folder,filesep);
home_folder = home_folder(1:slash_idx(end-1));
load([home_folder, 'locations.mat'],'-mat')
if ~exist(locations.scope)
    error(['Invalid mount location for images: "',locations.scope,...
        '" not found. Please load, update, and re-save "locations.mat"'])
end
[~,status] = urlread(locations.spreadsheet);
if ~status
    error(['Spreadsheet URL load unsuccessful. Please load, update, and re-save "locations.mat"'])
end

% Read in content from the "Scope Runs - Tracking Sets" spreadsheet
data = readScopeRuns(locations.spreadsheet, cell2mat(varargin));


% Display the sets we're running
for idx = 1:numel(data.dates)
    disp(['- ',data.dates{idx},'_',data.names{idx}])
end


% Cycle/measure sets
for idx = 1:numel(data.dates)
    
    % PARAMETERS
    parameters.SaveDirectory = [data.save_dir{idx},filesep,data.dates{idx},'_',data.names{idx}];
    load([parameters.SaveDirectory,filesep,'AllMeasurements.mat'])
    disp([parameters.SaveDirectory,', before: ', num2str(AllMeasurements.parameters.FramesPerHour),' frames per hr'])
    clear p;
    eval(data.modify{idx});
    if exist('p','var'); 
        AllMeasurements.parameters = combinestructures(p, AllMeasurements.parameters);
    else
        AllMeasurements.parameters.FramesPerHour = 12;
    end
    disp([parameters.SaveDirectory,', after: ', num2str(AllMeasurements.parameters.FramesPerHour),' frames per hr'])
    disp('- - - - -')
    save([parameters.SaveDirectory,filesep,'AllMeasurements.mat'],'AllMeasurements')
    clear AllMeasurements
end
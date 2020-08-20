function last_falltime = get_last_falltime(Trajectories,responder_index, varargin)
% get_last_falltime
%called 'computeDuration' in Ade's and Apeksha's code
% -----------------------------------------------------------------------
% INPUT:
%   REQUIRED:
%       Trajectories:
%   OPTIONAL:
%       Verbose: true or false
%       'Threshold'
p = inputParser;
addRequired(p, 'Trajectories', @isnumeric);
addRequired(p, 'responder_index');
addOptional(p, 'Verbose',false ,@islogical);
addParameter(p, 'Threshold',  0.2, @isnumeric);

parse(p, Trajectories, responder_index, varargin{:});
data = p.Results.Trajectories;
thresh = p.Results.Threshold;

convertFun = @(y) y/12; %frames to hour
neg = @(j) j*-1;
normFun = @(x) normalize(x, 'range');
last_falltime = zeros(size(data,1),1);

for i =1: size(data,1)
   %% all non-responders get NaN for last_falltime 
  %  if responder_index(i) == 0
   %     last_falltime (i) = nan;
    %else
        %%
        endpt =find(isfinite(data(i,:)),1,'last');
        X = fillmissing(data(i,1:endpt), 'linear')';
        t = transpose(1:endpt);
        f=fit(t,X , 'smoothingspline','SmoothingParam',0.2);
        x0 = feval(f, t);
        x1 = differentiate(f,t);

        %%  Find strongest decline rate
        z0 = normFun(x0);
        Z = normFun(zscore(neg(x1)));
        [~,zPkLoc] = findpeaks(Z,'MinPeakProminence',thresh, 'MinPeakWidth',3);
        if isempty(zPkLoc) % if too stringent, find peaks
            [~,zPkLoc]= findpeaks(Z);
        end
        if isempty(zPkLoc)
            last_falltime(i) = 0;
        else
            %% Find where trajectories taper off, i.e. derivative => 0
            pkIx = max(zPkLoc);
            % find where the signal crosses the base of the last peak prominence
            [~, maxDur]=getPeakBounds(z0,pkIx);
            last_falltime(i) = convertFun(maxDur);
        end
    %end
end
%% Plots
if p.Results.Verbose
    for i = randi([1 size(data, 1)], 1, 10)
    plot([0:1/12:(size(data, 2)-1)/12], data(i, :))
    hold on
    plot(last_falltime(i),1,'r*')
    keyboard
    clf
    end
end

function[leftBound, rightBound]= getPeakBounds(x,peakLoc)
%Returns the left boundary & right boundary of a peak
params = {'MinSeparation', 6, 'FlatSelection', 'first'};

leftBound=peakLoc - find(islocalmin(flip(x(1:peakLoc)), params{:}),1,'first')+1;
rightBound = find(islocalmin(x(peakLoc:end), params{:}),1, 'first') + peakLoc-1;
leftBound = max(leftBound, 1);
if isempty(rightBound)
    rightBound =  numel(x);
end
end
end
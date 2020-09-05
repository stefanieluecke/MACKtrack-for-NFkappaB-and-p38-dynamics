function [peaks, locs, HMW, prom, heights] = globalpeaks(vect, num_peaks)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% [peaks, locs] = globalpeaks(vect, num_peaks)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% GLOBALPEAKS seeks to find the "dominant" peaks in an input vector - output will be sorted
% from most-to-least dominant.
%
% INPUT:
% vect            input vector
% num_peaks       number of desired "dominant" peaks
%
% OUTPUT:
% peaks          peak values - are outputed as baseline subtracted
% locs           peak locations (index)
% heights        peak height (above nearest troughs)
% HMW            peak width
% prom           peak prominence (similar to height, except now includes
% height above nearest troughs where all neighboring peaks are shorter)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

% Find all peaks in vector; 1st global peak is maximum value overall
[all_peaks, all_locs, all_HMW, all_prom] = findpeaks(vect, 'WidthReference', 'halfheight');
all_peaks(((all_locs==1)) | (all_locs==length(vect)) ) = [];
all_HMW(((all_locs==1)) | (all_locs==length(vect)) ) = [];
all_prom(((all_locs==1)) | (all_locs==length(vect)) ) = [];
all_locs(((all_locs==1)) | (all_locs==length(vect)) ) = [];


getclosest = @(idx,vect) vect(find(abs(vect-idx)==min(abs(vect-idx)),1,'first'));
peaks = [];
locs  = [];
HMW = [];
prom = [];
heights = [];

while length(peaks) < num_peaks
    % Eliminate peaks that have been found already
    tmp_peaks = all_peaks;
    tmp_locs = all_locs;
    tmp_HMW = all_HMW;
    tmp_prom = all_prom;
    tmp_peaks(ismember(tmp_locs,locs)) = [];
    tmp_locs(ismember(tmp_locs,locs)) = [];
    tmp_HMW(ismember(tmp_locs,locs)) = [];
    tmp_prom(ismember(tmp_locs,locs)) = [];
    if isempty(tmp_peaks)
      break
    end
   
    % For each candidate, identify nearest peaks - maximize difference btw candidate and two nearest troughs.  
    diffs = zeros(size(tmp_peaks));
    loc_compare = [1 locs length(vect)];
    for i = 1:length(tmp_locs)
        tmp = loc_compare; tmp(tmp>=tmp_locs(i)) = inf;
        trough1 = min(vect(getclosest(tmp_locs(i),tmp):tmp_locs(i)));
        tmp = loc_compare; tmp(tmp<=tmp_locs(i)) = inf;
        trough2 = min(vect(tmp_locs(i):getclosest(tmp_locs(i),tmp)));
        diffs(i) = tmp_peaks(i) - max([trough1, trough2]);
    end
%todo understand prom and height difference properly

    peaks = [peaks, tmp_peaks(find(diffs==max(diffs),1,'first'))];
    locs = [locs, tmp_locs(find(diffs==max(diffs),1,'first'))];
    HMW = [HMW, tmp_HMW(find(diffs==max(diffs),1,'first'))];
    prom = [prom, tmp_prom(find(diffs==max(diffs),1,'first'))];
    heights = [heights, max(diffs)];
end
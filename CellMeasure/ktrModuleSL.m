function [CellMeasurements, ModuleDataOut] = ktrModuleSL(CellMeasurements,parameters, labels, AuxImages, ModuleData)

% NOT USED YET! TESTING (testing sth, maybe comparing to NFKB module, I
% don't remember...)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% KTRMODULESL  Attempt to make KTR Module more similar to NFkBdim module, % testing 

%              measure cytoplasmic:nuclear ratio in auxiliary (fluorescent) image
%
% CellMeasurements    structure with fields corresponding to cell measurements
%
% parameters          experiment data (total cells, total images, output directory)
% labels              Cell,Nuclear label matricies (labels.Cell and labels.Nucleus)
% AuxImages           images to measure
% ModuleData          extra information (current ModuleData.iter, etc.) used in measurement 
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% 

%IMAGE CORRECTION % 3 possibilities for Background correction
% Background correct, method 1: use flatfield image. 
if length(AuxImages)>1 % If more than one image in the fluorescence channel?? If more than one fluorescence channel?  
    if isequal(size(AuxImages{1}),size(parameters.Flatfield{1})) %check whether size of AuxImage and Flatfield image provided match %I think we provide a flatfield image in the Macktrack GUI
        AuxImages{1} = double(AuxImages{1}) - double(parameters.Flatfield{end}); %create double-precision arrays, substract last flatfield image from first AuxImage???
        ktr = flatfieldcorrect(AuxImages{1},double(parameters.Flatfield{1}));%nfkb is defined as flatfield corrected AuxImage pixels 
    else
        error(['Size mismatch between provided flatfield (#', num2str(img), ' and AuxImage'])
    end
else
% Flatfield correct, method 2: apply quadratic model to image
    if ~isfield(ModuleData,'X') %do not fully understand when this is used?
        ModuleData.X = backgroundcalculate(size(AuxImages{1}));
    end
    warning off MATLAB:nearlySingularMatrix
    pStar = (ModuleData.X'*ModuleData.X)\(ModuleData.X')*double(AuxImages{1}(:));
    warning on MATLAB:nearlySingularMatrix
    % Apply background correction
    ktr = reshape((double(AuxImages{1}(:) - ModuleData.X*pStar)),size(AuxImages{1}));
    ktr = ktr-min(ktr(:)); % Set minimum to zero
end

% Background correct (use unimodal model)
if ~isfield(ModuleData,'distr') %do not fully understand when this is used?
    [~, ModuleData.distr] = modebalance(ktr,1,ModuleData.BitDepth,'measure');
else
    ktr = modebalance(ktr,1,ModuleData.BitDepth,'correct',ModuleData.distr);
end


% MEASUREMENT
% On first call, initialize all new CellMeasurements fields 
 if ~isfield(CellMeasurements,['KTR_ratio',num2str(i)])
    CellMeasurements.(['ktr_ratio',num2str(i)]) =  nan(parameters.TotalCells,parameters.TotalImages,1);
    CellMeasurements.(['ktr_nuc',num2str(i)]) =  nan(parameters.TotalCells,parameters.TotalImages,1);
    CellMeasurements.(['ktr_cyto',num2str(i)]) = nan(parameters.TotalCells,parameters.TotalImages,1);
    
    
% Count cells for each frame, initialize bins
cells = unique(labels.Nucleus);
cells(cells==0) = [];

%%% This is present in NFkBdim, but not NFkB module --> needed here?
%[~, distr2] = modebalance(nfkb(labels.Cell==0),1,ModuleData.BitDepth,'measure'); %some form of background correction

% Cycle through each image and assign measurements
for i = 1:length(cells)
    nucleus = labels.Nucleus==cells(i);
    cytoplasm = (labels.Cell==cells(i)) &~nucleus;
    
     median_n = median(ktr(nucleus));
    %consider testing mode, mean, median for cytoplasmic fraction (see nfkbModule)???
    %add position correction using modebalance for cytoplasm and/or nucleus before moving on to create ratio?
     median_c = median(ktr(cytoplasm)); 
% Assign measurements
CellMeasurements.ktr_nuc(cells(i),ModuleData.iter) = median_n(1); %what is the 1 for?
CellMeasurements.ktr_cyto(cells(i),ModuleData.iter) = median_c(1); %test median and mean and mode 
CellMeasurements.ktr_ratio(cells(i),ModuleData.iter) = median_c(1)/median_c(1);
% Clear all variables, close all figure and files
clear
close all
fclose all;
% Add Binning folder and subfolders to search path
addpath(genpath('F:\Users\fmami\Documents\TUDelft\Thesis\LargeScaleCodes\Binning'))
addpath(genpath('F:\Users\fmami\Documents\TUDelft\Thesis\LargeScaleCodes\MitrottaThesis'))

%% Ask whether to load old .mat file for binning
answer = questdlg('Load old .mat file for binning?',...
    'Load .mat file','Yes','No','No');
switch answer
    case 'Yes'
        % Select .mat file to load
        [matFilename,matFilepath] = uigetfile('*.mat',...
            'Select one .mat file');
        load(strcat(matFilepath,matFilename))
        % Retrieve binning input parameters
        binningInputStruct = binnningInputParameters;
        
    case 'No'
        %% Retrieve binning input parameters
        binningInputStruct = binnningInputParameters;
        %% Set parameters for phase averaging
        % Path to folder with files containing a reference signal for phase
        % averaging
        referenceSignalFilesFolderPath = 'F:\Users\fmami\Documents\TUDelft\Thesis\ExperimentNovember\LabView\Wind tunnel test 2018_11_26\Data';
        % Index vector that maps the LabVIEW batches to the StbRun objects
        labview2StbRunIndexVector = [1:4,6:7,8:14,16];
%         labview2StbRunIndexVector = 1:2;
        % Nominal frequency of the studied phenomenon
        nominalGustFrequency = 1;  % [Hz]
        
        %% Select files to load
        % Select STB .dat files to be loaded
        stbDataPathArray = selectStbDatFiles(strcat(...
            binningInputStruct.projectFolderPath,...
            binningInputStruct.stbFolderName));
        % Select the files containing the reference signals to be loaded
        labviewDataPathArray = selectLabviewFiles(...
            referenceSignalFilesFolderPath);
        
        %% Load STB data and other files
        % Generate StbRun object
        stbRunVector = StbRun(stbDataPathArray);
        % Generate LabviewData object
        labviewBatchArray = LabviewData(labviewDataPathArray);
        
        %% Find reference initial time and cycle period for each STB run
        % Iterate through the different STB runs
        for i=length(stbRunVector):-1:1
            fprintf(['Parsing LabVIEW batches for identification of ',...
                'cycles and phases [%d/%d]\n'],length(stbRunVector)-i+1,...
                length(stbRunVector))
            % Calculate the period and the instant of the first phase of the
            % reference periodic signal from LabVIEW data
            [phase0ReferenceTime,referenceGustPeriod,~,pivTimeVector] =...
                labviewBatchArray(labview2StbRunIndexVector(i)...
                ).calculatePhaseAverageReference(nominalGustFrequency);
            % Store reference gust period
            referencePeriodVector(i) = referenceGustPeriod;
            % Calculate the time offset to apply to the STB data
            timeOffsetVector(i) = pivTimeVector(1)-phase0ReferenceTime;
        end
        
        %% Find binning input variables
        [binningInputArray,particleNondimensionalTimeVector] =...
            stbRunVector.getPhaseAverageBinningInput(...
            referencePeriodVector,timeOffsetVector);
        % Save the nondimensional time vector into binningInputStruct
        binningInputStruct.particleNondimensionalTimeVector=...
            particleNondimensionalTimeVector;
        % Set phase average flag to true
        binningInputStruct.phaseAverageFlag = true;
        
        %% Save workspace
        k = strfind(stbDataPathArray{1},filesep);
        saveFolderPath = stbDataPathArray{1}(1:k(end));
        filename = matlab.lang.makeValidName(['PreBinning',datestr(now,...
            'yyyymmdd-HHMM')]);
        save([saveFolderPath,filename,'.mat'],'stbRunVector',...
            'binningInputStruct','referencePeriodVector',...
            'timeOffsetVector','binningInputArray','saveFolderPath',...
            '-v7.3')
end

%% Perform binning
binnedData = binning(binningInputArray,binningInputStruct.binSize,...
    binningInputStruct);

%% Save data
% Save results in .mat file
fprintf('Saving binning results in .mat file\n')
binSizeString = sprintf('%.1fmm%.3fs',...
    binningInputStruct.binSize*1e3,...
    binningInputStruct.temporalBinSize);
filename = matlab.lang.makeValidName([...
    binningInputStruct.averagingMethod,'-BinSize',...
    binSizeString,'-BinnedOn',datestr(now,'yyyymmdd-HHMM')]);
save([saveFolderPath,filename,'.mat'],'binningInputArray',...
    'binningInputStruct','binnedData','-v7.3')
% Write binned data to tecplot file
fprintf('Writing binary file for tecplot\n')
binnedData.writeTecplotBinary(filename,saveFolderPath)

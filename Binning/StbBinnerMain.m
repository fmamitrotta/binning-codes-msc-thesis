% Clear all variables, close all figure and files
clear
close all
fclose all;
% Add Binning folder and subfolders to search path (substitute string with 
% your local path)
addpath(genpath(...
    'F:\Users\fmami\Documents\TUDelft\Thesis\LargeScaleCodes\Binning'))

%% Retrieve binning input parameters
binningInputStruct = binnningInputParameters;

%% Select files to load
% Select STB .dat files to be loaded
stbDataPathArray = selectStbDatFiles(strcat(...
    binningInputStruct.projectFolderPath,...
    binningInputStruct.stbFolderName));
% Select any other file here
%...........................

%% Load STB data and other files
% Generate StbRun object
stbRunArray = StbRun(stbDataPathArray);
% Load any other file here
%.........................

%% Additional code
% This part can be used to calculate the reference period and the time
% offset for phase average binning
%.................................

%% Perform binning
% Use timeAverageBinning or phaseAverageBinning according to what you want
% to do
[binningInputArray,binnedData] = stbRunArray.timeAverageBinning(...
    binningInputStruct);

%% Save data
filename = [stbDataPathArray{1}(1:end-4),'-',matlab.lang.makeValidName([...
    'BinnedOn',datestr(now,'yyyymmdd-HHMM')])];
save([filename,'.mat'],...
    'binningInputArray','binnedData',...
    'binningInputStruct','stbDataPathArray')

%% Write binned data to tecplot file
binnedData.writeTecplotBinary(filename)

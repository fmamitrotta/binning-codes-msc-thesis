function labviewDataPathArray = selectLabviewFiles(projectFolder)

%%  Compile list with runs to be binned %%

% Select .txt files to load
[labviewFileArray,labviewPath] = uigetfile('*.txt',...
    'Select One or More LabVIEW Files',...
    projectFolder,...
    'MultiSelect', 'on');
% Assemble cell array with complete path to the selected .dat files
if ischar(labviewFileArray)
    % If only one .dat file has been selected generate cell array
    labviewDataPathArray = {horzcat(labviewPath,labviewFileArray)};
else
    labviewDataPathArray = cellfun(@(x) horzcat(labviewPath,x),...
        labviewFileArray,'UniformOutput',false)';
end
end

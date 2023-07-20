function stbDataPathArray = selectStbDatFiles(projectFolder)

%%  Compile list with runs to be binned %%

% Select .dat files to load
[stbFileArray,stbPath] = uigetfile('*.dat',...
    'Select One or More STB Files',...
    projectFolder,...
    'MultiSelect', 'on');
% Assemble cell array with complete path to the selected .dat files
if ischar(stbFileArray)
    % If only one .dat file has been selected generate cell array
    stbDataPathArray = {horzcat(stbPath,stbFileArray)};
else
    stbDataPathArray = cellfun(@(x) horzcat(stbPath,x),stbFileArray,...
        'UniformOutput',false)';
end
end

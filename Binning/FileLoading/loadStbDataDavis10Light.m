function [particleXyzUvwArray,particleTimeArray] =...
    loadStbDataDavis10Light(stbFilePathArray)

% Define parameters for textscan
delimiter = {' ','='};
noHeaderLines = 3;
timeRowFormatSpec = '%*s%*s%*s%f%[^\n\r]';
dataRowFormatSpec = [repmat('%f',1,13),'%[^\n\r]'];
% Iterate through the input STB files
for i=length(stbFilePathArray):-1:1
    fprintf('Importing STB data from file [%d/%d]\n-> %s\n',...
        length(stbFilePathArray)-i+1,length(stbFilePathArray),...
        stbFilePathArray{i})
    % Open current STB file
    fileId = fopen(stbFilePathArray{i},'r');
    % Initialize counters
    j = 0;
    lastRow = 0;
    % Initialize acquisition arrays
    particleTimeArray{i,1} = NaN(1e7,1);
    rawXyzUvwArray{i,1} = NaN(1e7,6);
    % Read first time step and the related set of data
    partialTimeVector = textscan(fileId,timeRowFormatSpec,1,...
        'HeaderLines',noHeaderLines,'Delimiter',delimiter);
    partialDataArray = textscan(fileId,dataRowFormatSpec,...
        'HeaderLines',noHeaderLines,'Delimiter',delimiter);
    %     tmp2 = textscan(fileId,formatSpec,'HeaderLines',startRow,'Delimiter',delimiter);
    % Iterate until the end of file
    while (size(partialDataArray{1},1)>=1)
        % Increment time step counter
        j = j+1;
        % Determine the number of particles for current time step
        noParticlesRead = size(partialDataArray{1},1);
        % Store acquired time instant
        particleTimeArray{i}(lastRow+1:lastRow+noParticlesRead,:) =...
            partialTimeVector{1}(1)*ones(noParticlesRead,1);
        % Store acquired positional and velocity information
        rawXyzUvwArray{i}(lastRow+1:lastRow+noParticlesRead,:) = [...
            partialDataArray{1},partialDataArray{2},partialDataArray{3},...
            partialDataArray{5},partialDataArray{6},partialDataArray{7}];
        % Update the index of last row of fullDataArray
        lastRow = lastRow+noParticlesRead;
        % Read new time step and related set of data
        partialTimeVector = textscan(fileId,timeRowFormatSpec,1,...
            'HeaderLines',noHeaderLines-2,'Delimiter',delimiter);
        partialDataArray = textscan(fileId,dataRowFormatSpec,...
            'HeaderLines',noHeaderLines,'Delimiter',delimiter);
    end
    % Close current STB file
    fclose(fileId);
    % Clean acquisition arrays from nan
    particleTimeArray{i} = particleTimeArray{i}(~isnan(...
        particleTimeArray{i}));
    rawXyzUvwArray{i} = reshape(rawXyzUvwArray{i}(~isnan(...
        rawXyzUvwArray{i})),lastRow,6);
end
% Assembe output array with positional and velocity information
particleXyzUvwArray = cell2mat(rawXyzUvwArray);
% Output positional information in meters
particleXyzUvwArray(:,1:3) = particleXyzUvwArray(:,1:3)*1e-3;
end

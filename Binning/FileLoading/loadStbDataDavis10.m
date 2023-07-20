function [particleVector,trackVector,timeStepObjVector] =...
    loadStbDataDavis10(stbFilepath)

% Delimiter parameter for textscan
delimiter = {' ','='};
% Set number of lines to be skipped for the different calls of textscan
noStartHeaderLines = 2;
time2DataLineOffset = 3;
% Set format for time instants and rest of data
timeFormatSpec = '%*s%*s%*s%f%[^\n\r]';
dataArrayFormatSpec = [repmat('%f',1,13),'%[^\n\r]'];
% Initialize acquisition arrays
noParticlesPerTimeStepVector = NaN(1e7,1);
timeStepVector = NaN(1e7,1);
timeInstantVector = NaN(1e7,1);
fullDataArray  = NaN(1e7,13);
% Initialize counters
i = 0;
dataArrayLastRow = 0;
% Open STB file
fileId = fopen(stbFilepath,'r');
% Read first time step and the related set of data
partialTimeVector = textscan(fileId,timeFormatSpec,2,...
    'HeaderLines',noStartHeaderLines,'Delimiter',delimiter);
partialDataArray = textscan(fileId,dataArrayFormatSpec,...
    'HeaderLines',time2DataLineOffset,'Delimiter',delimiter);
% Iterate until the end of file
while (size(partialDataArray{1},1)>=1)
    % Increment time step counter
    i = i+1;
    % Determine the number of particles for current time step
    noParticlesRead = size(partialDataArray{1},1);
    noParticlesPerTimeStepVector(i) = noParticlesRead;
    % Update time vector and data array with info of current time step
    timeStepVector(i) = partialTimeVector{1}(1);
    timeInstantVector(i) = partialTimeVector{1}(2);
    fullDataArray(...
        dataArrayLastRow+1:dataArrayLastRow+noParticlesRead,:) =...
        [partialDataArray{1},partialDataArray{2},partialDataArray{3},...
        partialDataArray{4},partialDataArray{5},partialDataArray{6},...
        partialDataArray{7},partialDataArray{8},partialDataArray{9},...
        partialDataArray{10},partialDataArray{11},partialDataArray{12},...
        partialDataArray{13}];
    % Update the index of last row of fullDataArray
    dataArrayLastRow = dataArrayLastRow+noParticlesRead;
    % Read new time step and related set of data
    partialTimeVector = textscan(fileId,timeFormatSpec,2,...
        'Delimiter',delimiter);
    partialDataArray = textscan(fileId,dataArrayFormatSpec,...
        'HeaderLines',time2DataLineOffset,'Delimiter',delimiter);
end
% Close file
fclose(fileId);
% Clean acquisition arrays from nan
noParticlesPerTimeStepVector = noParticlesPerTimeStepVector(~isnan(...
    noParticlesPerTimeStepVector));
timeStepVector = timeStepVector(~isnan(timeStepVector));
timeInstantVector = timeInstantVector(~isnan(timeInstantVector));
fullDataArray = reshape(fullDataArray(~isnan(fullDataArray)),...
    dataArrayLastRow,13);
% Columns in dataArray represent the following
% 1 x[mm]
% 2 y[mm]
% 3 z[mm]
% 4 I
% 5 u[m/s]
% 6 v[m/s]
% 7 w[m/s]
% 8 |V|[m/s]
% 9 trackID
% 10 ax[m/s²]
% 11 ay[m/s²]
% 12 az[m/s²]
% 13 |a|[m/s²]
% Convert position units from mm to m
fullDataArray(:,1:3) = fullDataArray(:,1:3)*1e-3;
% Generate Particle object
particleVector = StbParticle(fullDataArray(:,1),...
    fullDataArray(:,2),...
    fullDataArray(:,3),...
    fullDataArray(:,4),...
    fullDataArray(:,5),...
    fullDataArray(:,6),...
    fullDataArray(:,7),...
    fullDataArray(:,8),...
    fullDataArray(:,10),...
    fullDataArray(:,11),...
    fullDataArray(:,12),...
    fullDataArray(:,13));
% Assemble Particle objects per time step
noParticlesPerTimeStepVector = [0;noParticlesPerTimeStepVector];
particlePerTimeStepArray = arrayfun(@(x) particleVector(...
    sum(noParticlesPerTimeStepVector(1:x))+1:...
    sum(noParticlesPerTimeStepVector(1:x))+...
    noParticlesPerTimeStepVector(x+1)),...
    1:length(noParticlesPerTimeStepVector)-1,'UniformOutput',false)';
% Generate TimeStep object
timeStepObjVector = StbTimeStep(timeStepVector,timeInstantVector,...
    particlePerTimeStepArray);
% Assemble Particle objects per track id
% trackIdVector includes the track id of all particles, as a consequence
% particles from the same track generate duplicates in the vector. Once the
% vector is obtained, it is necessary to map the indices of same track id
% to the correct track in order to assemble the Track object with the
% correct children Particle objects
trackIdVector = fullDataArray(:,9);
[particle2TrackIndexVector,uniqueTrackIdVector] =...
    findgroups(trackIdVector);
particlePerTrackArray = splitapply(@(x){particleVector(x)},...
    (1:length(trackIdVector))',particle2TrackIndexVector);
% Generate Track object
trackVector = StbTrack(num2cell(uniqueTrackIdVector),...
    particlePerTrackArray);
end

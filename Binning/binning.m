function binnedDataObject = binning(particleXyzUvwArray,binSize,...
    optionStruct)
%binning Perform binning on a STB dataset and return a BinnedStbData object
%   binnedDataObject = binning(particleXyzUvwArray,xLim,yLim,zLim,binSize,
%   optionStruct) carries out an ensemble average of the Lagrangian
%   particle tracks obtained by Shake-the-Box and return a BinnedStbData
%   object containing the eulerian description of the flow.
%   particleXyzUvwArray is a double array having the positional and
%   velocity information of the particles to bin in each row (columns are
%   x,y,z coordinates and u,v,w velocity components. xLim, yLim and zLim
%   are vectors with two elements each indicating the limits of the volume
%   to be binned. binSize is the dimension of the spatial bin and
%   optionStruct is a structure containing all kinds of optional input for
%   the ensemble average process (such as the overlap factor, the minimum
%   number of particles, the averaging method, options for phase average
%   and pressure calculation, etc.).
% Determin binning volume limits
if isfield(optionStruct,'domainCropVector')
    % If option to crop binning volume is specified in the
    % input structure, then set limits accordingly
    xLim = [optionStruct.domainCropVector(1),...
        optionStruct.domainCropVector(2)];
    yLim = [optionStruct.domainCropVector(3),...
        optionStruct.domainCropVector(4)];
    zLim = [optionStruct.domainCropVector(5),...
        optionStruct.domainCropVector(6)];
else
    % If no option to crop the binning volume is specified,
    % then set the limits considering the minimum and maximum
    % values of the particles' coordinates
    xLim = [min(particleXyzUvwArray(:,1)),max(particleXyzUvwArray(:,1))];
    yLim = [min(particleXyzUvwArray(:,2)),max(particleXyzUvwArray(:,2))];
    zLim = [min(particleXyzUvwArray(:,3)),max(particleXyzUvwArray(:,3))];
end
% Set the default overlap factor
overlapFactor = 0.75;
% If provided by the user, modify the overlap factor
if nargin>=3 && isfield(optionStruct,'overlapFactor')
    overlapFactor = optionStruct.overlapFactor;
end
% Calculate the location of bins' centroid according to the volume limits
% and the vector spacing (given by bin size and overlap factor). Take last
% centroids in each direction exceeding the volume limits by one vector
% spacing length in order to compensate for the next operation
vectorSpacing = binSize*(1-overlapFactor);
binCentroidXvector = xLim(1):vectorSpacing:xLim(2)+vectorSpacing;
binCentroidYvector = yLim(1):vectorSpacing:yLim(2)+vectorSpacing;
binCentroidZvector = zLim(1):vectorSpacing:zLim(2)+vectorSpacing;
% Shift the location of the bins' centroid such that the coordinates of the
% first centroid are a multiple of the spatial bin size (convenient for
% visualization
binCentroidXvector = binCentroidXvector-mod(binCentroidXvector(1),...
    vectorSpacing);
binCentroidYvector = binCentroidYvector-mod(binCentroidYvector(1),...
    vectorSpacing);
binCentroidZvector = binCentroidZvector-mod(binCentroidZvector(1),...
    vectorSpacing);
% Set default phase average flag
phaseAverageFlag = false;
% If provided by the user, modify the averaging method
if nargin>=3 && isfield(optionStruct,'phaseAverageFlag')
    phaseAverageFlag = optionStruct.phaseAverageFlag;
end
% If phase average binning, then find the particles that may be included in
% temporal bins spanning between the end and the start of the average cycle
if optionStruct.phaseAverageFlag
    % Find indices of particles located within half temporal bin from the
    % start of the average cycle. These particles have to be repeated at
    % the end of the average cycle so that temporal bins centred towards
    % the end of the cycle can include such particles
    appendLogicalVector =...
        optionStruct.particleNondimensionalTimeVector<=...
        optionStruct.temporalBinSize/2;
    % Find indices of particles located within half temporal bin from the
    % end of the average cycle. These particles have to be repeated at
    % the start of the average cycle so that temporal bins centred towards
    % the start of the cycle can include such particles
    prependLogicalVector = 1-...
        optionStruct.particleNondimensionalTimeVector<=...
        optionStruct.temporalBinSize/2;
    % Repeat positional and velocity information of identified particles
    particleXyzUvwArray = [particleXyzUvwArray(prependLogicalVector,:);...
        particleXyzUvwArray;particleXyzUvwArray(appendLogicalVector,:)];
    % Add corresponding nondimensional time
    optionStruct.particleNondimensionalTimeVector = [...
        optionStruct.particleNondimensionalTimeVector(...
        prependLogicalVector)-1;...
        optionStruct.particleNondimensionalTimeVector;...
        optionStruct.particleNondimensionalTimeVector(...
        appendLogicalVector)+1];
end
% Divide the volume of interest into a number of x slices
% Set the default number of slices
noSlices = 15;
% If provided by the user, modify the averaging method
if nargin>=3 && isfield(optionStruct,'noSlices')
    noSlices = optionStruct.noSlices;
end
% Partition the column containing the x coordinates into a matrix
% containing in each column nonoverlapping data segments. The last column
% is zero-padded to reach a coherent length with respect to other columns
binCentroidSliceXarray = buffer(binCentroidXvector,max([round(length(...
    binCentroidXvector)/noSlices),1]));
% Convert the matrix into a cell array with each cell containing a matrix
% column
binCentroidSliceXarray = num2cell(binCentroidSliceXarray,1);
% Remove the trailing zeros from the last cell of the array
binCentroidSliceXarray{end} = binCentroidSliceXarray{end}(1:find(...
    binCentroidSliceXarray{end},1,'last'));
% Find the index vector indicating the arrangement of the particles orederd
% by increasing x position
[~,I] = sort(particleXyzUvwArray(:,1));
% Assemble an array containing all positional and velocity information of
% the particles ordered by increasing x position
particleXyzUvwArray = particleXyzUvwArray(I,:);
% Set the search distance for particles' position form the centroid of the
% bins as half of the spatial bin size
searchDistance = binSize/2;
% Set the default minimum number of particles
minNoParticles = 5;
% If provided by the user, modify the averaging method
if nargin>=3 && isfield(optionStruct,'minNoParticles')
    minNoParticles = optionStruct.minNoParticles;
end
% Set the default averaging method
averagingMethod = 'adaptive';
% If provided by the user, modify the averaging method
if nargin>=3 && isfield(optionStruct,'averagingMethod')
    averagingMethod = optionStruct.averagingMethod;
end
% Distinguish between time average and phase average binning
if ~phaseAverageFlag
    %% Time average binning
    % Define function handle for binning according to the desired averaging
    % method
    fprintf(['Time average binning, selecting ensemble average method',...
        '\n-> %s\n'],averagingMethod)
    switch averagingMethod
        case 'tophat'
            binningFunction = @tophatAverage;
        case 'gaussian'
            binningFunction = @gaussianTimeAverage;
        case 'linear'
            binningFunction = @linearTimeAverage;
        case 'quadratic'
            binningFunction = @quadraticTimeAverage;
        case 'adaptive'
            binningFunction = @adaptiveTimeAverage;
        otherwise
            % If no valid averaging method is provided, issue a warning a
            % return control to invoking function
            warning(['Invalid averaging method, binning stopped.\nUse',...
                'one among ''tophat'', ''gaussian'', ''linear'', ',...
                '''quadratic'' and ''adaptive'''])
            return
    end
    % Look for the interval of particleXyzUvWArray indices for
    % which the particles are in the slice, using the binarySearch
    % function. binarySearch is a mex function that performs the binary
    % search algorithm to find "item(s)" (the values to be searched
    % for) among some pre-sorted "data" vector.
    % pos = binarySearch(data, items) where data is the pre-sorted
    % vector, items the desired values and pos is index of the first
    % instance of each item or the index of the closest item if the
    % items are not found
    for i=size(binCentroidSliceXarray,2):-1:1
        % First particle within the x interval of the current slice
        lowerIndex = binarySearch(particleXyzUvwArray(:,1),...
            binCentroidSliceXarray{i}(1)-searchDistance);
        % Last particle within the x interval of the current slice
        upperIndex = binarySearch(particleXyzUvwArray(:,1),...
            binCentroidSliceXarray{i}(end)+searchDistance);
        % Extract positional and velocity information of particles
        % belonging to current slice
        particleSliceXyzUvwArray{i} = particleXyzUvwArray(...
            lowerIndex:upperIndex,:);
    end
    % Perform binning of slices
    [outputVelocityMeanArray,outputVelocityStdArray,...
        outputNoFoundParticlesArray,outputUprimeUprimeArray,...
        outputVprimeVprimeArray,outputWprimeWprimeArray,...
        outputUprimeVprimeArray,outputUprimeWprimeArray,...
        outputVprimeWprimeArray] = spatialSliceBinning(...
        binCentroidSliceXarray,binCentroidYvector,binCentroidZvector,...
        particleSliceXyzUvwArray,minNoParticles,binningFunction,binSize);
    % At the end of the iteration through the slices of the binned volume
    % the output cell arrays are transformed into double arrays collapsing
    % the data from different slices together
    outputVelocityMeanArray = cell2mat(outputVelocityMeanArray);
    outputVelocityStdArray = cell2mat(outputVelocityStdArray);
    outputNoFoundParticlesArray = cell2mat(outputNoFoundParticlesArray);
    outputUprimeUprimeArray = cell2mat(outputUprimeUprimeArray);
    outputVprimeVprimeArray = cell2mat(outputVprimeVprimeArray);
    outputWprimeWprimeArray = cell2mat(outputWprimeWprimeArray);
    outputUprimeVprimeArray = cell2mat(outputUprimeVprimeArray);
    outputUprimeWprimeArray = cell2mat(outputUprimeWprimeArray);
    outputVprimeWprimeArray = cell2mat(outputVprimeWprimeArray);
    % Generate the output BinnedStdData object
    binnedDataObject = BinnedStbData(binCentroidXvector,...
        binCentroidYvector,binCentroidZvector,outputVelocityMeanArray,...
        outputVelocityStdArray,outputNoFoundParticlesArray,...
        outputUprimeUprimeArray,outputVprimeVprimeArray,...
        outputWprimeWprimeArray,outputUprimeVprimeArray,...
        outputUprimeWprimeArray,outputVprimeWprimeArray);
else
    %% Phase average binning
    % Check the use of temporal binning
    useTemporalBinning = true;
    if nargin>=3 && isfield(optionStruct,'useTemporalBinning')
        useTemporalBinning = optionStruct.useTemporalBinning;
    end
    % Define function handle for binning according to the desired averaging
    % method and distinguishing between the use of temporal binning or not
    fprintf(['Phase average binning, selecting ensemble average method',...
        '\n-> %s\n'],averagingMethod)
    if useTemporalBinning
        sliceBinningFunction = @spatioTemporalSliceBinning;
        switch averagingMethod
            case 'tophat'
                binningFunction = @tophatAverage;
            case 'gaussian'
                binningFunction = @gaussianPhaseAverage;
            case 'linear'
                binningFunction = @linearPhaseAverage;
            case 'quadratic'
                binningFunction = @quadraticPhaseAverage;
            case 'adaptive'
                binningFunction = @adaptivePhaseAverage;
            otherwise
                % If no valid averaging method is provided, issue a warning a
                % return control to invoking function
                warning(['Invalid averaging method, binning stopped.\nUse ',...
                    'one among ''tophat'', ''gaussian'', ''linear'', ',...
                    '''quadratic'' and ''adaptive'''])
                return
        end
    else
        sliceBinningFunction = @spatialSliceBinning;
        switch averagingMethod
            case 'tophat'
                binningFunction = @tophatAverage;
            case 'gaussian'
                binningFunction = @gaussianTimeAverage;
            case 'linear'
                binningFunction = @linearTimeAverage;
            case 'quadratic'
                binningFunction = @quadraticTimeAverage;
            case 'adaptive'
                binningFunction = @adaptiveTimeAverage;
            otherwise
                % If no valid averaging method is provided, issue a warning a
                % return control to invoking function
                warning(['Invalid averaging method, binning stopped.\nUse ',...
                    'one among ''tophat'', ''gaussian'', ''linear'', ',...
                    '''quadratic'' and ''adaptive'''])
                return
        end
    end
    % Order the nondimensional time vector provided in optionStruct in the
    % same say as done for positional and velocity information of the
    % particles
    particleNondimensionalTimeVector =...
        optionStruct.particleNondimensionalTimeVector(I);
    % Define the centroid of the temporal bins
    binNondimensionalTimeVector = linspace(0,1,optionStruct.noPhases+1)';
    % Store temporal bin size in variable in order to avoid the use of
    % optionStruct in the parfor
    temporalBinSize = optionStruct.temporalBinSize;
    % Iterate through the phases
    for t = length(binNondimensionalTimeVector):-1:1
        fprintf('Phase %d/%d\n',length(binNondimensionalTimeVector)-t+1,...
            length(binNondimensionalTimeVector))
        % Select particles within the temporal bin corresponding to the
        % current phase
        phaseSelectionLogical = particleNondimensionalTimeVector>=...
            binNondimensionalTimeVector(t)-temporalBinSize/2 &...
            particleNondimensionalTimeVector<=...
            binNondimensionalTimeVector(t)+temporalBinSize/2;
        particlePhaseNondimensionalTimeVector =...
            particleNondimensionalTimeVector(phaseSelectionLogical);
        particlePhaseXyzUvwArray =...
            particleXyzUvwArray(phaseSelectionLogical,:);
        % Look for the interval of particlePhaseXyzUvwArray indices for
        % which the particles are in the slice, using the binarySearch
        % function. binarySearch is a mex function that performs the binary
        % search algorithm to find "item(s)" (the values to be searched
        % for) among some pre-sorted "data" vector.
        % pos = binarySearch(data, items) where data is the pre-sorted
        % vector, items the desired values and pos is index of the first
        % instance of each item or the index of the closest item if the
        % items are not found
        for i=size(binCentroidSliceXarray,2):-1:1
            % First particle within the x interval of the current slice
            lowerIndex = binarySearch(particlePhaseXyzUvwArray(:,1),...
                binCentroidSliceXarray{i}(1)-searchDistance);
            % Last particle within the x interval of the current slice
            upperIndex = binarySearch(particlePhaseXyzUvwArray(:,1),...
                binCentroidSliceXarray{i}(end)+searchDistance);
            % Extract positional and velocity information of particles
            % belonging to current slice
            particleSliceXyzUvwArray{i} = particlePhaseXyzUvwArray(...
                lowerIndex:upperIndex,:);
            % Extract temporal information of particles belonging to
            % current slice
            particleSliceNondimensionalTimeVector{i} =...
                particlePhaseNondimensionalTimeVector(lowerIndex:...
                upperIndex,:);
        end
        % Perform binning of slices
        [outputVelocityMeanArray,outputVelocityStdArray,...
            outputNoFoundParticlesArray,outputUprimeUprimeArray,...
            outputVprimeVprimeArray,outputWprimeWprimeArray,...
            outputUprimeVprimeArray,outputUprimeWprimeArray,...
            outputVprimeWprimeArray] = sliceBinningFunction(...
            binCentroidSliceXarray,binCentroidYvector,...
            binCentroidZvector,particleSliceXyzUvwArray,minNoParticles,...
            binningFunction,binSize,binNondimensionalTimeVector(t),...
            particleSliceNondimensionalTimeVector,temporalBinSize);
        % At the end of the iteration through the slices of the binned volume
        % the output cell arrays are transformed into double arrays collapsing
        % the data from different slices together
        outputPhaseVelocityMeanArray{t} =...
            cell2mat(outputVelocityMeanArray);
        outputPhaseVelocityStdArray{t} = cell2mat(outputVelocityStdArray);
        outputPhaseNoFoundParticlesArray{t} =...
            cell2mat(outputNoFoundParticlesArray);
        outputPhaseUprimeUprimeArray{t} =...
            cell2mat(outputUprimeUprimeArray);
        outputPhaseVprimeVprimeArray{t} =...
            cell2mat(outputVprimeVprimeArray);
        outputPhaseWprimeWprimeArray{t} =...
            cell2mat(outputWprimeWprimeArray);
        outputPhaseUprimeVprimeArray{t} =...
            cell2mat(outputUprimeVprimeArray);
        outputPhaseUprimeWprimeArray{t} =...
            cell2mat(outputUprimeWprimeArray);
        outputPhaseVprimeWprimeArray{t} =...
            cell2mat(outputVprimeWprimeArray);
    end
    % Generate a vector of BinnedStdData objects representing the different
    % phases of the average cycle
    fprintf('Generating BinnedStbData object\n')
    binnedDataObject = BinnedStbData(binCentroidXvector,...
        binCentroidYvector,binCentroidZvector,...
        outputPhaseVelocityMeanArray,outputPhaseVelocityStdArray,...
        outputPhaseNoFoundParticlesArray,outputPhaseUprimeUprimeArray,...
        outputPhaseVprimeVprimeArray,outputPhaseWprimeWprimeArray,...
        outputPhaseUprimeVprimeArray,outputPhaseUprimeWprimeArray,...
        outputPhaseVprimeWprimeArray,binNondimensionalTimeVector);
end
end

%% Spatial binning of the volume slices
function [outputVelocityMeanArray,outputVelocityStdArray,...
    outputNoFoundParticlesArray,outputUprimeUprimeArray,...
    outputVprimeVprimeArray,outputWprimeWprimeArray,...
    outputUprimeVprimeArray,outputUprimeWprimeArray,...
    outputVprimeWprimeArray] = spatialSliceBinning(...
    binCentroidSliceXarray,binCentroidYvector,binCentroidZvector,...
    particleSliceXyzUvwArray,minNoParticles,binningFunction,binSize,...
    varargin)
% Iterate through the slices of the volume
% parfor i = 1:size(binCentroidSliceXarray,2)
for i = size(binCentroidSliceXarray,2):-1:1
    fprintf('Binning slice %d/%d\n',i,size(binCentroidSliceXarray,2))
    % Generate the matrix containing the 3D grid coordinates of the
    % current slice of the volume
    [xMeshSliceArray,yMeshSliceArray,zMeshSliceArray] = meshgrid(...
        binCentroidSliceXarray{i},binCentroidYvector,...
        binCentroidZvector);
    % Generate the matrix containing the 3D grid coordinates of the
    % current slice of the volume in columns |x|y|z|
    binCentroidSliceArray = [xMeshSliceArray(:),...
        yMeshSliceArray(:),zMeshSliceArray(:)];
    % rangesearch(X,Y,r) finds all the X points that are within distance r
    % of the Y points. The rows of X and Y correspond to observations, and
    % the columns correspond to variables. In this case we are finding all
    % the particleSliceXyzUvwArray locations that are within a distance of
    % half the linear bin size from the 3D grid centroids of the current
    % volume slice. indexArray is an my-by-1 cell array, where my is the
    % number of rows in binCentroidSliceArray. The cell indexArray{j}
    % contains the vector of indices of particleSliceXyzUvwArray points
    % whose distances from binCentroidSliceArray(j,:) are not greater than
    % searchDistance
    indexArray = rangesearch(particleSliceXyzUvwArray{i}(:,1:3),...
        binCentroidSliceArray,binSize/2,'distance',...
        'euclidean');
    % Initialize result arrays for current slice with NaN matrices.
    % This is done to output NaN variables for all centroids where no
    % particle or not enough particles are found to carry out the
    % ensemble average
    binnedVelocityMeanArray = NaN(length(indexArray),3);
    binnedVelocityStdArray = NaN(length(indexArray),3);
    noFoundParticlesArray = NaN(length(indexArray),1);
    uprimeUprimeArray = NaN(length(indexArray),1);
    uprimeVprimeArray = NaN(length(indexArray),1);
    uprimeWprimeArray = NaN(length(indexArray),1);
    vprimeVprimeArray = NaN(length(indexArray),1);
    vprimeWprimeArray = NaN(length(indexArray),1);
    wprimeWprimeArray = NaN(length(indexArray),1);
    % Iterate over centroids of current volume slice having at
    % least the minimum number of particles
    noParticlesInBinVector = cellfun(@(x) length(x),indexArray)';
    fprintf('-> %d bins with at least %d particles found\n',...
        length(find(noParticlesInBinVector>=minNoParticles)),...
        minNoParticles)
    for j = find(noParticlesInBinVector>=minNoParticles)
        % Extract positional and velocity information of particles
        % within the current spatial bin
        particleBinXyzArray = particleSliceXyzUvwArray{i}(...
            indexArray{j},1:3);
        particleBinUvwArray = particleSliceXyzUvwArray{i}(...
            indexArray{j},4:6);
        % Carry out binning ony for bins which coordinates are within the
        % min max interval of the particles' coordinates (this is done in
        % order to avoid extrapolation when polynomial fitting is used)
        if xMeshSliceArray(j)>=min(particleBinXyzArray(:,1)) &&...
                xMeshSliceArray(j)<=max(particleBinXyzArray(:,1)) &&...
                yMeshSliceArray(j)>=min(particleBinXyzArray(:,2)) &&...
                yMeshSliceArray(j)<=max(particleBinXyzArray(:,2)) &&...
                zMeshSliceArray(j)>=min(particleBinXyzArray(:,3)) &&...
                zMeshSliceArray(j)<=max(particleBinXyzArray(:,3))
            % Velocity filter: remove particles where all velocity
            % components are further than 3 time the standard deviation
            % from the mean
            velocityMeanArray = mean(particleBinUvwArray);
            velocityStdArray = std(particleBinUvwArray);
            filteringVector = sum(abs(particleBinUvwArray-...
                ones(size(particleBinUvwArray,1),1)*...
                velocityMeanArray)>3*ones(size(...
                particleBinUvwArray,1),1)*velocityStdArray,2);
            particleBinFilteredXyzArray =...
                particleBinXyzArray(filteringVector<3,:);
            particleBinFilteredUvwArray =...
                particleBinUvwArray(filteringVector<3,:);
            % Call binning function and return mean velocity, standard
            % deviation and velocity fluctuations
            [binVelocityMeanVector,binVelocityStdVector,...
                velocityFluctuationArray] = binningFunction(...
                particleBinFilteredXyzArray,particleBinFilteredUvwArray,...
                minNoParticles,xMeshSliceArray(j),yMeshSliceArray(j),...
                zMeshSliceArray(j),binSize);
            % Store results into array
            binnedVelocityMeanArray(j,:) = binVelocityMeanVector;
            binnedVelocityStdArray(j,:) = binVelocityStdVector;
            % Number of found particles corresponds to the number of
            % rows of particleOverallBinFilteredXyzArray
            noFoundParticlesArray(j,1) =...
                size(particleBinFilteredXyzArray,1);
            % Calculate Reynolds stresses from velocity
            % fluctuations
            uprimeUprimeArray(j,1) = mean(...
                velocityFluctuationArray(:,1).*...
                velocityFluctuationArray(:,1));
            uprimeVprimeArray(j,1) = mean(...
                velocityFluctuationArray(:,1).*...
                velocityFluctuationArray(:,2));
            uprimeWprimeArray(j,1) = mean(...
                velocityFluctuationArray(:,1).*...
                velocityFluctuationArray(:,3));
            vprimeVprimeArray(j,1) = mean(...
                velocityFluctuationArray(:,2)...
                .*velocityFluctuationArray(:,2));
            vprimeWprimeArray(j,1) = mean(...
                velocityFluctuationArray(:,2).*...
                velocityFluctuationArray(:,3));
            wprimeWprimeArray(j,1) = mean(...
                velocityFluctuationArray(:,3).*...
                velocityFluctuationArray(:,3));
        end
    end
    % Reshape variables obtained from current slice and store it in the
    % output cell array. Variables are reshaped in order to have the
    % same dimension of the xMeshSliceArray, becoming in this way 3D
    % double arrays (the mean and std velocity arrays become 4D double
    % arrays because of the three velocity components)
    outputVelocityMeanArray{i} = reshape(binnedVelocityMeanArray,...
        [size(xMeshSliceArray),3]);
    outputVelocityStdArray{i} = reshape(binnedVelocityStdArray,...
        [size(xMeshSliceArray),3]);
    outputNoFoundParticlesArray{i} = reshape(noFoundParticlesArray,...
        size(xMeshSliceArray));
    outputUprimeUprimeArray{i} = reshape(uprimeUprimeArray,size(...
        xMeshSliceArray));
    outputVprimeVprimeArray{i} = reshape(vprimeVprimeArray,size(...
        xMeshSliceArray));
    outputWprimeWprimeArray{i} = reshape(wprimeWprimeArray,size(...
        xMeshSliceArray));
    outputUprimeVprimeArray{i} = reshape(uprimeVprimeArray,size(...
        xMeshSliceArray));
    outputUprimeWprimeArray{i} = reshape(uprimeWprimeArray,size(...
        xMeshSliceArray));
    outputVprimeWprimeArray{i} = reshape(vprimeWprimeArray,size(...
        xMeshSliceArray));
    fprintf('-> slice %d/%d binned\n',i,size(binCentroidSliceXarray,2))
end
end

%% Spatio-temporal binning of the volume slices
function [outputVelocityMeanArray,outputVelocityStdArray,...
    outputNoFoundParticlesArray,outputUprimeUprimeArray,...
    outputVprimeVprimeArray,outputWprimeWprimeArray,...
    outputUprimeVprimeArray,outputUprimeWprimeArray,...
    outputVprimeWprimeArray] = spatioTemporalSliceBinning(...
    binCentroidSliceXarray,binCentroidYvector,binCentroidZvector,...
    particleSliceXyzUvwArray,minNoParticles,binningFunction,...
    spatialBinSize,binNondimensionalTime,...
    particleSliceNondimensionalTimeVector,temporalBinSize)
% Iterate through the slices of the volume
parfor i = 1:size(binCentroidSliceXarray,2)
    % for i = size(binCentroidSliceXarray,2):-1:1
    fprintf('Binning slice %d/%d (phase t/T=%.2f)\n',i,...
        size(binCentroidSliceXarray,2),binNondimensionalTime)
    % Generate the matrix containing the 3D grid coordinates of the
    % current slice of the volume
    [xMeshSliceArray,yMeshSliceArray,zMeshSliceArray] = meshgrid(...
        binCentroidSliceXarray{i},binCentroidYvector,...
        binCentroidZvector);
    % Generate the matrix containing the 3D grid coordinates of the
    % current slice of the volume in columns |x|y|z|
    binCentroidSliceArray = [xMeshSliceArray(:),...
        yMeshSliceArray(:),zMeshSliceArray(:)];
    % rangesearch(X,Y,r) finds all the X points that are within distance r
    % of the Y points. The rows of X and Y correspond to observations, and
    % the columns correspond to variables. In this case we are finding all
    % the particleSliceXyzUvwArray locations that are within a distance of
    % half the linear bin size from the 3D grid centroids of the current
    % volume slice. indexArray is an my-by-1 cell array, where my is the
    % number of rows in binCentroidSliceArray. The cell indexArray{j}
    % contains the vector of indices of particleSliceXyzUvwArray points
    % whose distances from binCentroidSliceArray(j,:) are not greater than
    % searchDistance
    indexArray = rangesearch(particleSliceXyzUvwArray{i}(:,1:3),...
        binCentroidSliceArray,spatialBinSize/2,'distance','euclidean');
    % Initialize result arrays for current slice with NaN matrices.
    % This is done to output NaN variables for all centroids where no
    % particle or not enough particles are found to carry out the
    % ensemble average
    binnedVelocityMeanArray = NaN(length(indexArray),3);
    binnedVelocityStdArray = NaN(length(indexArray),3);
    noFoundParticlesArray = NaN(length(indexArray),1);
    uprimeUprimeArray = NaN(length(indexArray),1);
    uprimeVprimeArray = NaN(length(indexArray),1);
    uprimeWprimeArray = NaN(length(indexArray),1);
    vprimeVprimeArray = NaN(length(indexArray),1);
    vprimeWprimeArray = NaN(length(indexArray),1);
    wprimeWprimeArray = NaN(length(indexArray),1);
    % Iterate over centroids of current volume slice having at
    % least the minimum number of particles
    noParticlesInBinVector = cellfun(@(x) length(x),indexArray)';
    fprintf('-> %d bins with at least %d particles found\n',...
        length(find(noParticlesInBinVector>=minNoParticles)),...
        minNoParticles)
    for j = find(noParticlesInBinVector>=minNoParticles)
        % Extract positional and velocity information of particles
        % within the current spatial bin
        particleBinXyzArray = particleSliceXyzUvwArray{i}(...
            indexArray{j},1:3);
        particleBinUvwArray = particleSliceXyzUvwArray{i}(...
            indexArray{j},4:6);
        % Extract temporal information of particles within the
        % current spatial bin
        particleBinNondimensionalTimeVector =...
            particleSliceNondimensionalTimeVector{i}(...
            indexArray{j});
        % Carry out binning ony for bins which coordinates are within the
        % min max interval of the particles' spatial and temporal 
        % coordinates (this is done in order to avoid extrapolation when 
        % polynomial fitting is used)
        if xMeshSliceArray(j)>=min(particleBinXyzArray(:,1)) &&...
                xMeshSliceArray(j)<=max(particleBinXyzArray(:,1)) &&...
                yMeshSliceArray(j)>=min(particleBinXyzArray(:,2)) &&...
                yMeshSliceArray(j)<=max(particleBinXyzArray(:,2)) &&...
                zMeshSliceArray(j)>=min(particleBinXyzArray(:,3)) &&...
                zMeshSliceArray(j)<=max(particleBinXyzArray(:,3)) &&...
                binNondimensionalTime>=min(...
                particleBinNondimensionalTimeVector) &&...
                binNondimensionalTime<=max(...
                particleBinNondimensionalTimeVector)
            % Velocity filter: remove particles where all velocity
            % components are further than 3 time the standard deviation
            % from the mean
            velocityMeanArray = mean(particleBinUvwArray);
            velocityStdArray = std(particleBinUvwArray);
            filteringVector = sum(abs(particleBinUvwArray-...
                ones(size(particleBinUvwArray,1),1)*...
                velocityMeanArray)>3*ones(size(...
                particleBinUvwArray,1),1)*velocityStdArray,2);
            particleBinFilteredXyzArray =...
                particleBinXyzArray(filteringVector<3,:);
            particleBinFilteredUvwArray =...
                particleBinUvwArray(filteringVector<3,:);
            particleBinFilteredNondimensionalTimeVector =...
                particleBinNondimensionalTimeVector(...
                filteringVector<3,:);
            % Call binning function and return mean velocity, standard
            % deviation and velocity fluctuations
            [binVelocityMeanVector,binVelocityStdVector,...
                velocityFluctuationArray] = binningFunction(...
                particleBinFilteredXyzArray,...
                particleBinFilteredUvwArray,minNoParticles,...
                particleBinFilteredNondimensionalTimeVector,...
                xMeshSliceArray(j),yMeshSliceArray(j),...
                zMeshSliceArray(j),binNondimensionalTime,...
                spatialBinSize,temporalBinSize);
            % Store results into array
            binnedVelocityMeanArray(j,:) = binVelocityMeanVector;
            binnedVelocityStdArray(j,:) = binVelocityStdVector;
            % Number of found particles corresponds to the number of
            % rows of particleOverallBinFilteredXyzArray
            noFoundParticlesArray(j,1) =...
                size(particleBinFilteredXyzArray,1);
            % Calculate Reynolds stresses from velocity
            % fluctuations
            uprimeUprimeArray(j,1) = mean(...
                velocityFluctuationArray(:,1).*...
                velocityFluctuationArray(:,1));
            uprimeVprimeArray(j,1) = mean(...
                velocityFluctuationArray(:,1).*...
                velocityFluctuationArray(:,2));
            uprimeWprimeArray(j,1) = mean(...
                velocityFluctuationArray(:,1).*...
                velocityFluctuationArray(:,3));
            vprimeVprimeArray(j,1) = mean(...
                velocityFluctuationArray(:,2)...
                .*velocityFluctuationArray(:,2));
            vprimeWprimeArray(j,1) = mean(...
                velocityFluctuationArray(:,2).*...
                velocityFluctuationArray(:,3));
            wprimeWprimeArray(j,1) = mean(...
                velocityFluctuationArray(:,3).*...
                velocityFluctuationArray(:,3));
        end
    end
    % Reshape variables obtained from current slice and store it in the
    % output cell array. Variables are reshaped in order to have the
    % same dimension of the xMeshSliceArray, becoming in this way 3D
    % double arrays (the mean and std velocity arrays become 4D double
    % arrays because of the three velocity components)
    outputVelocityMeanArray{i} = reshape(binnedVelocityMeanArray,...
        [size(xMeshSliceArray),3]);
    outputVelocityStdArray{i} = reshape(binnedVelocityStdArray,...
        [size(xMeshSliceArray),3]);
    outputNoFoundParticlesArray{i} = reshape(noFoundParticlesArray,...
        size(xMeshSliceArray));
    outputUprimeUprimeArray{i} = reshape(uprimeUprimeArray,size(...
        xMeshSliceArray));
    outputVprimeVprimeArray{i} = reshape(vprimeVprimeArray,size(...
        xMeshSliceArray));
    outputWprimeWprimeArray{i} = reshape(wprimeWprimeArray,size(...
        xMeshSliceArray));
    outputUprimeVprimeArray{i} = reshape(uprimeVprimeArray,size(...
        xMeshSliceArray));
    outputUprimeWprimeArray{i} = reshape(uprimeWprimeArray,size(...
        xMeshSliceArray));
    outputVprimeWprimeArray{i} = reshape(vprimeWprimeArray,size(...
        xMeshSliceArray));
    fprintf('-> slice %d/%d binned\n',i,size(binCentroidSliceXarray,2))
end
end

%% Time or phase averaged binning based on top-hat filter
function [binVelocityMeanVector,binVelocityStdVector,...
    velocityFluctuationArray] = tophatAverage(...
    particleBinFilteredXyzArray,particleBinFilteredUvwArray,...
    minNoParticles,~,~,~,~,~,~,~)
%tophatAverage Carry out an ensemble average in one bin by means of a
%top-hat filter (corresponding to a simple average)
%   Detailed explanation goes here
if size(particleBinFilteredXyzArray,1)>=minNoParticles
    % If number of particle is larger or equal than the minimum allowed,
    % then average with top-hat filter (simple ensemble average)
    % Calculate mean velocity
    binVelocityMeanVector = mean(particleBinFilteredUvwArray);
    % Calculate standard deviation with gaussian weights
    binVelocityStdVector = std(particleBinFilteredUvwArray);
    % Calculate velocity fluctuations
    velocityFluctuationArray = particleBinFilteredUvwArray-...
        binVelocityMeanVector;
else
    % If number of particles if smaller than the minimum
    % allowed, set all output variables to NaN
    binVelocityMeanVector = NaN(1,3);
    binVelocityStdVector = NaN(1,3);
    velocityFluctuationArray = NaN(size(particleBinFilteredXyzArray));
end
end

%% Time average binning based on Gaussian filter
function [binVelocityMeanVector,binVelocityStdVector,...
    velocityFluctuationArray] = gaussianTimeAverage(...
    particleBinFilteredXyzArray,particleBinFilteredUvwArray,...
    minNoParticles,binXcoordinate,binYcoordinate,binZcoordinate,...
    spatialBinSize)
%gaussianTimeAverage Carry out an time ensemble average in one bin by means
%of a Gaussian filter.
%   Detailed explanation goes here
if size(particleBinFilteredXyzArray,1)>=minNoParticles
    % If number of particle is larger or equal than the
    % minimum allowed, then average with gaussian filter
    [binVelocityMeanVector,binVelocityStdVector,...
        velocityFluctuationArray] = spatialGaussianFilter(...
        particleBinFilteredXyzArray,particleBinFilteredUvwArray,...
        binXcoordinate,binYcoordinate,binZcoordinate,spatialBinSize);
else
    % If number of particles if smaller than the minimum
    % allowed, set all output variables to NaN
    binVelocityMeanVector = NaN(1,3);
    binVelocityStdVector = NaN(1,3);
    velocityFluctuationArray = NaN(size(particleBinFilteredXyzArray));
end
end

%% Time average binning based on linear polynomial fit
function [binVelocityMeanVector,binVelocityStdVector,...
    velocityFluctuationArray] = linearTimeAverage(...
    particleBinFilteredXyzArray,particleBinFilteredUvwArray,...
    minNoParticles,binXcoordinate,binYcoordinate,binZcoordinate,~)
%linearTimeAverage Carry out a time ensemble average in one bin by means of
%a linear polynomial fit.
%   Detailed explanation goes here
if size(particleBinFilteredXyzArray,1)>=max([minNoParticles,8])
    % If number of particles is larger or equal than the
    % minimum allowed and double the amount of coefficients
    % to be found, then average with a linear polynomial
    % fit
    [binVelocityMeanVector,binVelocityStdVector,...
        velocityFluctuationArray] = spatialLinearFit(...
        particleBinFilteredXyzArray,particleBinFilteredUvwArray,...
        binXcoordinate,binYcoordinate,binZcoordinate);
else
    % If number of particles if smaller than the minimum
    % allowed, set all output variables to NaN
    binVelocityMeanVector = NaN(1,3);
    binVelocityStdVector = NaN(1,3);
    velocityFluctuationArray = NaN(size(particleBinFilteredXyzArray));
end
end

%% Time average binning based on quadratic polynomial fit
function [binVelocityMeanVector,binVelocityStdVector,...
    velocityFluctuationArray] = quadraticTimeAverage(...
    particleBinFilteredXyzArray,particleBinFilteredUvwArray,...
    minNoParticles,binXcoordinate,binYcoordinate,binZcoordinate,~)
%quadraticTimeAverage Carry out a time ensemble average in one bin by means
%of a quadratic polynomial fit.
%   Detailed explanation goes here
if size(particleBinFilteredXyzArray,1)>=max([minNoParticles,20])
    % If number of particles is larger or equal than the
    % minimum allowed and double the amount of coefficients
    % to be found, then average with a linear polynomial
    % fit
    [binVelocityMeanVector,binVelocityStdVector,...
        velocityFluctuationArray] = spatialQuadraticFit(...
        particleBinFilteredXyzArray,particleBinFilteredUvwArray,...
        binXcoordinate,binYcoordinate,binZcoordinate);
else
    % If number of particles if smaller than the minimum
    % allowed, set all output variables to NaN
    binVelocityMeanVector = NaN(1,3);
    binVelocityStdVector = NaN(1,3);
    velocityFluctuationArray = NaN(size(particleBinFilteredXyzArray));
end
end

%% Adaptive time average binning
function [binVelocityMeanVector,binVelocityStdVector,...
    velocityFluctuationArray] = adaptiveTimeAverage(...
    particleBinFilteredXyzArray,particleBinFilteredUvwArray,...
    minNoParticles,binXcoordinate,binYcoordinate,binZcoordinate,...
    spatialBinSize)
%adaptiveTimeAverage Carry out a time ensemble average in one bin adapting
%the method to the number of particles available.
%   Detailed explanation goes here
if size(particleBinFilteredXyzArray,1)>=max([minNoParticles,20])
    % If number of particles is larger or equal than the
    % minimum allowed and double the amount of coefficients
    % to be found, then average with a linear polynomial
    % fit
    [binVelocityMeanVector,binVelocityStdVector,...
        velocityFluctuationArray] = spatialQuadraticFit(...
        particleBinFilteredXyzArray,particleBinFilteredUvwArray,...
        binXcoordinate,binYcoordinate,binZcoordinate);
    %
elseif size(particleBinFilteredXyzArray,1)>=max(...
        [minNoParticles,8])
    % If number of particles is larger or equal than the
    % minimum allowed and double the amount of coefficients
    % to be found, then average with a linear polynomial
    % fit
    [binVelocityMeanVector,binVelocityStdVector,...
        velocityFluctuationArray] = spatialLinearFit(...
        particleBinFilteredXyzArray,particleBinFilteredUvwArray,...
        binXcoordinate,binYcoordinate,binZcoordinate);
    %
elseif size(particleBinFilteredXyzArray,1)>=...
        minNoParticles
    % If number of particle is larger or equal than the
    % minimum allowed, then average with gaussian filter
    [binVelocityMeanVector,binVelocityStdVector,...
        velocityFluctuationArray] = spatialGaussianFilter(...
        particleBinFilteredXyzArray,particleBinFilteredUvwArray,...
        binXcoordinate,binYcoordinate,binZcoordinate,spatialBinSize);
else
    % If number of particles if smaller than the minimum
    % allowed, set all output variables to NaN
    binVelocityMeanVector = NaN(1,3);
    binVelocityStdVector = NaN(1,3);
    velocityFluctuationArray = NaN(size(particleBinFilteredXyzArray));
end
end

%% Phase average binning based on Gaussian filter
function [binVelocityMeanVector,binVelocityStdVector,...
    velocityFluctuationArray] = gaussianPhaseAverage(...
    particleBinFilteredXyzArray,particleBinFilteredUvwArray,...
    minNoParticles,particleBinFilteredNondimensionalTimeVector,...
    binXcoordinate,binYcoordinate,binZcoordinate,binNondimensionalTime,...
    spatialBinSize,temporalBinSize)
%gaussianPhaseAverage Carry out an phase ensemble average in one
%spatio-temporal bin by means of a spatio-temporal Gaussian filter.
%   Detailed explanation goes here
if size(particleBinFilteredXyzArray,1)>=minNoParticles
    % If number of particle is larger or equal than the
    % minimum allowed, then average with gaussian filter
    [binVelocityMeanVector,binVelocityStdVector,...
        velocityFluctuationArray] = spatioTemporalGaussianFilter(...
        particleBinFilteredXyzArray,particleBinFilteredUvwArray,...
        particleBinFilteredNondimensionalTimeVector,binXcoordinate,...
        binYcoordinate,binZcoordinate,binNondimensionalTime,...
        spatialBinSize,temporalBinSize);
else
    % If number of particles if smaller than the minimum
    % allowed, set all output variables to NaN
    binVelocityMeanVector = NaN(1,3);
    binVelocityStdVector = NaN(1,3);
    velocityFluctuationArray = NaN(size(particleBinFilteredXyzArray));
end
end

%% Phase average binning based on linear polynomial fit
function [binVelocityMeanVector,binVelocityStdVector,...
    velocityFluctuationArray] = linearPhaseAverage(...
    particleBinFilteredXyzArray,particleBinFilteredUvwArray,...
    minNoParticles,particleBinFilteredNondimensionalTimeVector,...
    binXcoordinate,binYcoordinate,binZcoordinate,binNondimensionalTime,...
    ~,~)
%linearPhaseAverage Carry out an phase ensemble average in one
%spatio-temporal bin by means of a spatio-temporal linear polynomial fit.
%   Detailed explanation goes here
if size(particleBinFilteredXyzArray,1)>=max([minNoParticles,10])
    % If number of particles is larger or equal than the
    % minimum allowed and double the amount of coefficients
    % to be found, then average with a linear polynomial
    % fit
    [binVelocityMeanVector,binVelocityStdVector,...
        velocityFluctuationArray] = spatioTemporalLinearFit(...
        particleBinFilteredXyzArray,particleBinFilteredUvwArray,...
        particleBinFilteredNondimensionalTimeVector,binXcoordinate,...
        binYcoordinate,binZcoordinate,binNondimensionalTime);
else
    % If number of particles if smaller than the minimum
    % allowed, set all output variables to NaN
    binVelocityMeanVector = NaN(1,3);
    binVelocityStdVector = NaN(1,3);
    velocityFluctuationArray = NaN(size(particleBinFilteredXyzArray));
end
end

%% Phase average binning based on quadratic polynomial fit
function [binVelocityMeanVector,binVelocityStdVector,...
    velocityFluctuationArray] = quadraticPhaseAverage(...
    particleBinFilteredXyzArray,particleBinFilteredUvwArray,...
    minNoParticles,particleBinFilteredNondimensionalTimeVector,...
    binXcoordinate,binYcoordinate,binZcoordinate,binNondimensionalTime,...
    ~,~)
%quadraticPhaseAverage Carry out an phase ensemble average in one
%spatio-temporal bin by means of a spatio-temporal quadratic polynomial
%fit.
%   Detailed explanation goes here
if size(particleBinFilteredXyzArray,1)>=max([minNoParticles,30])
    % If number of particles is larger or equal than the
    % minimum allowed and double the amount of coefficients
    % to be found, then average with a linear polynomial
    % fit
    [binVelocityMeanVector,binVelocityStdVector,...
        velocityFluctuationArray] = spatioTemporalQuadraticFit(...
        particleBinFilteredXyzArray,particleBinFilteredUvwArray,...
        particleBinFilteredNondimensionalTimeVector,binXcoordinate,...
        binYcoordinate,binZcoordinate,binNondimensionalTime);
else
    % If number of particles if smaller than the minimum
    % allowed, set all output variables to NaN
    binVelocityMeanVector = NaN(1,3);
    binVelocityStdVector = NaN(1,3);
    velocityFluctuationArray = NaN(size(particleBinFilteredXyzArray));
end
end

%% Adaptive phase average binning
function [binVelocityMeanVector,binVelocityStdVector,...
    velocityFluctuationArray] = adaptivePhaseAverage(...
    particleBinFilteredXyzArray,particleBinFilteredUvwArray,...
    minNoParticles,particleBinFilteredNondimensionalTimeVector,...
    binXcoordinate,binYcoordinate,binZcoordinate,binNondimensionalTime,...
    spatialBinSize,temporalBinSize)
%adaptivePhaseAverage  Carry out a time ensemble average in one
%spatio-temporal bin adapting the method to the number of particles
%available.
%   Detailed explanation goes here
if size(particleBinFilteredXyzArray,1)>=max([minNoParticles,30])
    % If number of particles is larger or equal than the
    % minimum allowed and double the amount of coefficients
    % to be found, then average with a quadratic polynomial
    % fit
    [binVelocityMeanVector,binVelocityStdVector,...
        velocityFluctuationArray] = spatioTemporalQuadraticFit(...
        particleBinFilteredXyzArray,particleBinFilteredUvwArray,...
        particleBinFilteredNondimensionalTimeVector,binXcoordinate,...
        binYcoordinate,binZcoordinate,binNondimensionalTime);
    %
elseif size(particleBinFilteredXyzArray,1)>=max([minNoParticles,10])
    % If number of particles is larger or equal than the
    % minimum allowed and double the amount of coefficients
    % to be found, then average with a linear polynomial
    % fit
    [binVelocityMeanVector,binVelocityStdVector,...
        velocityFluctuationArray] = spatioTemporalLinearFit(...
        particleBinFilteredXyzArray,particleBinFilteredUvwArray,...
        particleBinFilteredNondimensionalTimeVector,binXcoordinate,...
        binYcoordinate,binZcoordinate,binNondimensionalTime);
    %
elseif size(particleBinFilteredXyzArray,1)>=minNoParticles
    % If number of particle is larger or equal than the
    % minimum allowed, then average with gaussian filter
    [binVelocityMeanVector,binVelocityStdVector,...
        velocityFluctuationArray] = spatioTemporalGaussianFilter(...
        particleBinFilteredXyzArray,particleBinFilteredUvwArray,...
        particleBinFilteredNondimensionalTimeVector,binXcoordinate,...
        binYcoordinate,binZcoordinate,binNondimensionalTime,...
        spatialBinSize,temporalBinSize);
else
    % If number of particles if smaller than the minimum
    % allowed, set all output variables to NaN
    binVelocityMeanVector = NaN(1,3);
    binVelocityStdVector = NaN(1,3);
    velocityFluctuationArray = NaN(size(particleBinFilteredXyzArray));
end
end

%% Spatial gaussian filter
function [binVelocityMeanVector,binVelocityStdVector,...
    velocityFluctuationArray] = spatialGaussianFilter(particleXyzArray,...
    particleUvwArray,binXcoordinate,binYcoordinate,binZcoordinate,...
    spatialBinSize)
%spatialGaussianFilter  Gaussian filter for ensemble averaging inside a
%spatial bin.
%   [binVelocityMeanVector,binVelocityStdVector,velocityFluctuationArray] =
%   spatialGaussianFilter(distanceArray,spatialBinNo,spatialBinSize,
%   particleOverallBinFilteredUvwArray) returns the mean velocity vector,
%   the related standard deviation and the constant mean based velocity
%   fluctuations resulting from the ensemble averaging with a gaussian
%   filter inside a spatial bin. Both mean and standard deviation are
%   obtained using gaussian weights based on the particle distances from
%   the bin centroid.
% Particle distances for spatial component of the weights
distance4GaussianWeightsVector = vecnorm(particleXyzArray-ones(size(...
    particleXyzArray,1),1)*[binXcoordinate,binYcoordinate,...
    binZcoordinate],2,2);
% Calculate gaussian weights
gaussianWeightVector = exp(-(distance4GaussianWeightsVector/...
    spatialBinSize).^2);
% Calculate mean velocity with gaussian weights
binVelocityMeanVector = sum(particleUvwArray.*(gaussianWeightVector*...
    ones(1,3)),1)./(sum(gaussianWeightVector)*ones(1,3));
% Calculate standard deviation with gaussian weights
binVelocityStdVector = std(particleUvwArray,gaussianWeightVector);
% Calculate velocity fluctuations
velocityFluctuationArray = particleUvwArray-binVelocityMeanVector;
end

%% Spatial linear polynomial fit
function [binVelocityMeanVector,binVelocityStdVector,...
    velocityFluctuationArray] = spatialLinearFit(particleXyzArray,...
    particleUvwArray,binXcoordinate,binYcoordinate,binZcoordinate)
%spatialLinearFit  Linear polynomial fit for ensemble averaging inside a
%spatial bin.
%   [binVelocityMeanVector,binVelocityStdVector,velocityFluctuationArray] =
%   spatialLinearFit(particleOverallBinFilteredXyzArray,xMeshSliceArray,
%   yMeshSliceArray,zMeshSliceArray,spatialBinNo,
%   particleOverallBinFilteredUvwArray) returns the mean velocity vector,
%   the related standard deviation and the varying mean based velocity
%   fluctuations resulting from the ensemble averaging with a linear
%   polynomial fit inside a spatial bin. The mean is based on a least
%   square regression of the velocity components inside the spatial bin,
%   which aims at predicting how the mean velocity varies inside the bin,
%   thus providing better estimate of velocity mean, standard deviation and
%   fluctuation.
% Calculate distances of particles from bin centroid
particleXyzDistanceArray = particleXyzArray-ones(size(particleXyzArray,...
    1),1)*[binXcoordinate,binYcoordinate,binZcoordinate];
% Assemble linear system of equations
A  = [ones(size(particleXyzDistanceArray,1),1),...
    particleXyzDistanceArray];
% Iterate over the three velocity components
for k = 3:-1:1
    % Find least-squares solution of the system of
    % equations
    fitCoefficientVector = A\particleUvwArray(:,k);
    % Retrieve mean velocity at bin centroid
    binVelocityMeanVector(k) = fitCoefficientVector(1);
    % Calculate the velocities inside the bin at the
    % positions and nondimensional time instants of the
    % fitted particles
    particleFittedUvwArray(:,k) = fitCoefficientVector(1)+...
        fitCoefficientVector(2)*particleXyzDistanceArray(:,1)+...
        fitCoefficientVector(3)*particleXyzDistanceArray(:,2)+...
        fitCoefficientVector(4)*particleXyzDistanceArray(:,3);
    % Calculate velocity standard deviation taking into
    % account the velocity variation inside the bin
    binVelocityStdVector(k) = std(particleUvwArray(:,k)-...
        particleFittedUvwArray(:,k));
end
% Calculate velocity fluctuations
velocityFluctuationArray = particleUvwArray-particleFittedUvwArray;
end

%% Spatial quadratic polynomial fit
function [binVelocityMeanArray,binVelocityStdArray,...
    velocityFluctuationArray] = spatialQuadraticFit(particleXyzArray,...
    particleUvwArray,binXcoordinate,binYcoordinate,binZcoordinate)
%spatialQuadraticFit  Quadratic polynomial fit for ensemble averaging
%inside a spatial bin.
%   [binVelocityMeanVector,binVelocityStdVector,velocityFluctuationArray] =
%   spatialLinearFit(particleOverallBinFilteredXyzArray,xMeshSliceArray,
%   yMeshSliceArray,zMeshSliceArray,spatialBinNo,
%   particleOverallBinFilteredUvwArray) returns the mean velocity vector,
%   the related standard deviation and the varying mean based velocity
%   fluctuations resulting from the ensemble averaging with a quadratic
%   polynomial fit inside a spatial bin. The mean is based on a least
%   square regression of the velocity components inside the spatial bin,
%   which aims at predicting how the mean velocity varies inside the bin,
%   thus providing better estimate of velocity mean, standard deviation and
%   fluctuation.
% Calculate distances of particles from bin centroid
particleXyzDistanceArray = particleXyzArray-ones(size(particleXyzArray,...
    1),1)*[binXcoordinate,binYcoordinate,binZcoordinate];
% Assemble linear system of equations
A  = [ones(size(particleXyzDistanceArray,1),1),...
    particleXyzDistanceArray,...
    particleXyzDistanceArray(:,1).^2,...
    particleXyzDistanceArray(:,1).*particleXyzDistanceArray(:,2),...
    particleXyzDistanceArray(:,2).^2,...
    particleXyzDistanceArray(:,1).*particleXyzDistanceArray(:,3),...
    particleXyzDistanceArray(:,2).*particleXyzDistanceArray(:,3),...
    particleXyzDistanceArray(:,3).^2];
% Iterate over the three velocity components
for k = 3:-1:1
    % Find least-squares solution of the system of
    % equations
    fitCoefficientVector = A\particleUvwArray(:,k);
    % Retrieve mean velocity at bin centroid
    binVelocityMeanArray(k) = fitCoefficientVector(1);
    % Calculate the velocities inside the bin at the
    % positions and nondimensional time instants of the
    % fitted particles
    particleFittedUvwArray(:,k) =...
        fitCoefficientVector(1)+...
        fitCoefficientVector(2)*particleXyzDistanceArray(:,1)+...
        fitCoefficientVector(3)*particleXyzDistanceArray(:,2)+...
        fitCoefficientVector(4)*particleXyzDistanceArray(:,3)+...
        fitCoefficientVector(5)*particleXyzDistanceArray(:,1).^2+...
        fitCoefficientVector(6)*particleXyzDistanceArray(:,1).*...
        particleXyzDistanceArray(:,2)+...
        fitCoefficientVector(7)*particleXyzDistanceArray(:,2).^2+...
        fitCoefficientVector(8)*particleXyzDistanceArray(:,1).*...
        particleXyzDistanceArray(:,3)+...
        fitCoefficientVector(9)*particleXyzDistanceArray(:,2).*...
        particleXyzDistanceArray(:,3)+...
        fitCoefficientVector(10)*particleXyzDistanceArray(:,3).^2;
    % Calculate velocity standard deviation taking into
    % account the velocity variation inside the bin
    binVelocityStdArray(k) = std(particleUvwArray(:,k)-...
        particleFittedUvwArray(:,k));
end
% Calculate velocity fluctuations
velocityFluctuationArray = particleUvwArray-particleFittedUvwArray;
end

%% Spatio-temporal gaussian filter
function [binVelocityMeanVector,binVelocityStdVector,...
    velocityFluctuationArray] = spatioTemporalGaussianFilter(...
    particleXyzArray,particleUvwArray,...
    particleNondimensionalTimeVector,binXcoordinate,...
    binYcoordinate,binZcoordinate,binNondimensionalTime,...
    spatialBinSize,temporalBinSize)
%spatioTemporalGaussianFilter  Gaussian filter for ensemble averaging
%inside a spatio-temporal bin.
%   [binVelocityMeanVector,binVelocityStdVector,velocityFluctuationArray] =
%   spatioTemporalGaussianFilter(distanceArray,spatialBinNo,timeIndexArray,
%   temporalBinNo,deltaTimeArray,spatialBinSize,temporalBinSize,
%   particleOverallBinFilteredUvwArray) returns the mean velocity vector,
%   the related standard deviation and the constant mean based velocity
%   fluctuations resulting from the ensemble averaging with a gaussian
%   filter inside a spatio-temporal bin. Both mean and standard deviation
%   are obtained using gaussian weights based on the particle distances
%   from the bin centroid and on the time delta from the instant of the
%   temporal bin.
% Particle distances for spatial component of the weights
distance4GaussianWeightsVector = vecnorm(particleXyzArray-ones(size(...
    particleXyzArray,1),1)*[binXcoordinate,binYcoordinate,...
    binZcoordinate],2,2);
% Particle time differences for temporal component of the weights
deltaTime4GaussianWeightsVector = particleNondimensionalTimeVector-...
    binNondimensionalTime;
% Calculate gaussian weights
gaussianWeightVector = exp(-((distance4GaussianWeightsVector/...
    spatialBinSize).^2+(deltaTime4GaussianWeightsVector/...
    temporalBinSize).^2));
% Calculate mean velocity with gaussian weights
binVelocityMeanVector = sum(particleUvwArray.*(gaussianWeightVector*...
    ones(1,3)),1)./(sum(gaussianWeightVector)*ones(1,3));
% Calculate standard deviation with gaussian weights
binVelocityStdVector = std(particleUvwArray,gaussianWeightVector);
% Calculate velocity fluctuations
velocityFluctuationArray = particleUvwArray-binVelocityMeanVector;
end

%% Spatio-temporal linear polynomial fit
function [binVelocityMeanVector,binVelocityStdVector,...
    velocityFluctuationArray] = spatioTemporalLinearFit(...
    particleXyzArray,particleUvwArray,...
    particleNondimensionalTimeVector,binXcoordinate,...
    binYcoordinate,binZcoordinate,binNondimensionalTime)
%spatioTemporalLinearFit  Linear polynomial fit for ensemble averaging
%inside a spatio-temporal bin.
%   [binVelocityMeanVector,binVelocityStdVector,velocityFluctuationArray] =
%   spatioTemporalLinearFit(particleXyzArray,particleUvwArray,
%   particleNondimensionalTimeVector,binXcoordinate,binYcoordinate,
%   binZcoordinate,binNondimensionalTime,spatialBinSize,temporalBinSize)
%   returns the mean velocity vector, the related standard deviation and
%   the varying mean based velocity fluctuations resulting from the
%   ensemble averaging with a linear polynomial fit inside a
%   spatio-temporal bin. The mean is based on a least square regression of
%   the velocity components inside the spatio-temporal bin, which aims at
%   predicting how the mean velocity varies in time and space inside the
%   bin, thus providing better estimate of velocity mean, standard
%   deviation and fluctuation.
% Calculate distances of particles from bin centroid
particleXyzDistanceArray = particleXyzArray-ones(size(particleXyzArray,...
    1),1)*[binXcoordinate,binYcoordinate,binZcoordinate];
% Calculate time difference of particles from the centre of the temporal
% bin
particleTimeDistanceVector = particleNondimensionalTimeVector-...
    binNondimensionalTime;
% Assemble linear system of equations
A  = [ones(size(particleXyzDistanceArray,1),1),...
    particleXyzDistanceArray,...
    particleTimeDistanceVector];
% Iterate over the three velocity components
for k = 3:-1:1
    % Find least-squares solution of the system of
    % equations
    fitCoefficientVector = A\particleUvwArray(:,k);
    % Retrieve mean velocity at bin centroid
    binVelocityMeanVector(k) = fitCoefficientVector(1);
    % Calculate the velocities inside the bin at the
    % positions and nondimensional time instants of the
    % fitted particles
    particleFittedUvwArray(:,k) = fitCoefficientVector(1)+...
        fitCoefficientVector(2)*particleXyzDistanceArray(:,1)+...
        fitCoefficientVector(3)*particleXyzDistanceArray(:,2)+...
        fitCoefficientVector(4)*particleXyzDistanceArray(:,3)+...
        fitCoefficientVector(5)*...
        particleTimeDistanceVector;
    % Calculate velocity standard deviation taking into
    % account the velocity variation inside the bin
    binVelocityStdVector(k) = std(particleUvwArray(:,k)-...
        particleFittedUvwArray(:,k));
end
% Calculate velocity fluctuations
velocityFluctuationArray = particleUvwArray-particleFittedUvwArray;
end

%% Spatio-temporal quadratic polynomial fit
function [binVelocityMeanArray,binVelocityStdArray,...
    velocityFluctuationArray] = spatioTemporalQuadraticFit(...
    particleXyzArray,particleUvwArray,...
    particleNondimensionalTimeVector,binXcoordinate,...
    binYcoordinate,binZcoordinate,binNondimensionalTime)
%spatioTemporalQuadraticFit  Quadratic polynomial fit for ensemble
%averaging inside a spatio-temporal bin.
%   [binVelocityMeanVector,binVelocityStdVector,velocityFluctuationArray] =
%   spatioTemporalQuadraticFit(particleXyzArray,particleUvwArray,
%   particleNondimensionalTimeVector,binXcoordinate,binYcoordinate,
%   binZcoordinate,binNondimensionalTime) returns the mean velocity vector,
%   the related standard deviation and the varying mean based velocity
%   fluctuations resulting from the ensemble averaging with a quadratic
%   polynomial fit inside a spatio-temporal bin. The mean is based on a
%   least square regression of the velocity components inside the
%   spatio-temporal bin, which aims at predicting how the mean velocity
%   varies in time and space inside the bin, thus providing better estimate
%   of velocity mean, standard deviation and fluctuation.
% Calculate distances of particles from bin centroid
particleXyzDistanceArray = particleXyzArray-ones(size(particleXyzArray,...
    1),1)*[binXcoordinate,binYcoordinate,binZcoordinate];
% Calculate time difference of particles from the centre of the temporal
% bin
particleTimeDistanceVector = particleNondimensionalTimeVector-...
    binNondimensionalTime;
% Assemble linear system of equations
A  = [ones(size(particleXyzDistanceArray,1),1),...
    particleXyzDistanceArray,...
    particleXyzDistanceArray(:,1).^2,...
    particleXyzDistanceArray(:,1).*particleXyzDistanceArray(:,2),...
    particleXyzDistanceArray(:,2).^2,...
    particleXyzDistanceArray(:,1).*particleXyzDistanceArray(:,3),...
    particleXyzDistanceArray(:,2).*particleXyzDistanceArray(:,3),...
    particleXyzDistanceArray(:,3).^2,...
    particleTimeDistanceVector,...
    particleTimeDistanceVector.*particleXyzDistanceArray(:,1),...
    particleTimeDistanceVector.*particleXyzDistanceArray(:,2),...
    particleTimeDistanceVector.*particleXyzDistanceArray(:,3),...
    particleTimeDistanceVector.^2];
% Iterate over the three velocity components
for k = 3:-1:1
    % Find least-squares solution of the system of
    % equations
    fitCoefficientVector = A\particleUvwArray(:,k);
    % Retrieve mean velocity at bin centroid
    binVelocityMeanArray(k) = fitCoefficientVector(1);
    % Calculate the velocities inside the bin at the
    % positions and nondimensional time instants of the
    % fitted particles
    particleFittedUvwArray(:,k) =...
        fitCoefficientVector(1)+...
        fitCoefficientVector(2)*particleXyzDistanceArray(:,1)+...
        fitCoefficientVector(3)*particleXyzDistanceArray(:,2)+...
        fitCoefficientVector(4)*particleXyzDistanceArray(:,3)+...
        fitCoefficientVector(5)*particleXyzDistanceArray(:,1).^2+...
        fitCoefficientVector(6)*particleXyzDistanceArray(:,1).*...
        particleXyzDistanceArray(:,2)+...
        fitCoefficientVector(7)*particleXyzDistanceArray(:,2).^2+...
        fitCoefficientVector(8)*particleXyzDistanceArray(:,1).*...
        particleXyzDistanceArray(:,3)+...
        fitCoefficientVector(9)*particleXyzDistanceArray(:,2).*...
        particleXyzDistanceArray(:,3)+...
        fitCoefficientVector(10)*particleXyzDistanceArray(:,3).^2+...
        fitCoefficientVector(11)*particleTimeDistanceVector+...
        fitCoefficientVector(12)*particleTimeDistanceVector.*...
        particleXyzDistanceArray(:,1)+...
        fitCoefficientVector(13)*particleTimeDistanceVector.*...
        particleXyzDistanceArray(:,2)+...
        fitCoefficientVector(14)*particleTimeDistanceVector.*...
        particleXyzDistanceArray(:,3)+...
        fitCoefficientVector(15)*particleTimeDistanceVector.^2;
    % Calculate velocity standard deviation taking into
    % account the velocity variation inside the bin
    binVelocityStdArray(k) = std(particleUvwArray(:,k)-...
        particleFittedUvwArray(:,k));
end
% Calculate velocity fluctuations
velocityFluctuationArray = particleUvwArray-particleFittedUvwArray;
end

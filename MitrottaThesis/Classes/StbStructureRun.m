classdef StbStructureRun < StbRun
    %StbStructureRun Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Dependent,SetAccess=private)
        MarkerVector    % vector of StructuralMarker objects
    end
    
    properties (Dependent)
        SearchRadius        % parameter for the union of tracks and markers
        FitTimeInterval     % parameter for the estimation of the markers' trajectory
        NoTimeStepsFilter   % parameter used to filter out markers
    end
    
    properties (Hidden,SetAccess=private)
        RootMarkerVector    % vector of unfiltered markers
        MarkerVectorProxy   % vector of filtered markers
        SearchRadiusProxy
        FitTimeIntervalProxy
        NoTimeStepsFilterProxy
    end
    
    methods
        %% Constructor
        function obj = StbStructureRun(stbFilepathArray,searchRadius,...
                fitTimeInterval,noTimeStepsFilter)
            %StbStructureRun Construct an instance of this class
            %   obj = StbStructureRun(stbFilepathArray,searchRadius,...
            %   fitTimeInterval,noTimeStepsFilter)
            % Set an empty path to input file is no input argument is given
            if nargin==0
                stbFilepathArray = {};
            end
            % Call superclass constructor
            obj = obj@StbRun(stbFilepathArray);
            % If 2 or more input arguments are given assign class specific
            % properties
            if nargin>=2
                % Assign search radius
                tempCellArray = num2cell(ones(size(obj))*searchRadius);
                [obj.SearchRadius] = tempCellArray{:};
                % Assign fit time interval
                tempCellArray = num2cell(ones(size(obj))*fitTimeInterval);
                [obj.FitTimeInterval] = tempCellArray{:};
                % If given as input variable assign the number of time
                % steps to filter markers
                if nargin==4
                    tempCellArray = num2cell(ones(size(obj))*...
                        noTimeStepsFilter);
                    [obj.NoTimeStepsFilter] = tempCellArray{:};
                end
            end
        end
        
        %% MarkerVector get method
        function markerVector = get.MarkerVector(obj)
            markerVector = obj.MarkerVectorProxy;
        end
        
        %% MarkerSearchRadius set and get method
        function set.SearchRadius(obj,searchRadius)
            % Set hidden property
            obj.SearchRadiusProxy = searchRadius;
            if ~isempty(obj.FitTimeInterval)
                % Generate the unfiltered vector of markers
                obj.RootMarkerVector = StructuralMarker(...
                    obj.TimeStepVector,obj.FitTimeInterval,...
                    searchRadius);
                % Filter out the markers with less than the prescribed
                % number of tracks
                if isempty(obj.NoTimeStepsFilter)
                    obj.MarkerVectorProxy = obj.RootMarkerVector;
                else
                    obj.MarkerVectorProxy =...
                        obj.RootMarkerVector.filterMarkersByNoTimeSteps(...
                        obj.NoTimeStepsFilter);
                end
            end
        end
        function searchRadius = get.SearchRadius(obj)
            searchRadius = obj.SearchRadiusProxy;
        end
        
        %% Time4MarkerFit set and get method
        function set.FitTimeInterval(obj,fitTimeInterval)
            % Set hidden property
            obj.FitTimeIntervalProxy = fitTimeInterval;
            if ~isempty(obj.SearchRadius)
                % Generate the unfiltered vector of markers
                obj.RootMarkerVector = StructuralMarker(...
                    obj.TimeStepVector,fitTimeInterval,obj.SearchRadius);
                % Filter out the markers with less than the prescribed
                % number of tracks
                if isempty(obj.NoTimeStepsFilter)
                    obj.MarkerVectorProxy = obj.RootMarkerVector;
                else
                    obj.MarkerVectorProxy =...
                        obj.RootMarkerVector.filterMarkersByNoTimeSteps(...
                        obj.NoTimeStepsFilter);
                end
            end
        end
        function fitTimeInterval = get.FitTimeInterval(obj)
            fitTimeInterval = obj.FitTimeIntervalProxy;
        end
        
        %% NoTimeStepsFilter set and get method
        function set.NoTimeStepsFilter(obj,noTimeStepsFilter)
            % Set hidden property
            obj.NoTimeStepsFilterProxy = noTimeStepsFilter;
            if ~isempty(obj.RootMarkerVector) && ~isempty(obj.SearchRadius)
                % Generate the filtered vector of markers
                obj.MarkerVectorProxy =...
                    obj.RootMarkerVector.filterMarkersByNoTimeSteps(...
                    noTimeStepsFilter);
            end
        end
        function noTimeStepsFilter = get.NoTimeStepsFilter(obj)
            noTimeStepsFilter = obj.NoTimeStepsFilterProxy;
        end
        
        %% Filter markers by number method
        function filterMarkersByNo(obj,noMarkers)
            %filterMarkersByNo Find the markers with largest number of time
            %steps
            %   filterMarkersByNo(obj,noMarkers) applies the number of time
            %   steps filter to find the noMarkers with the largest number
            %   of time steps
            % Iterate through the elements of the StbStructureRun vector
            % that have more markers than the parameter specified
            for i=find(noMarkers<arrayfun(@(x)...
                    length(x.RootMarkerVector),obj))'
                % Find number of time steps of all markers of current
                % StbStructureRun object
                noTimeStepsVector = arrayfun(@(x) length(x.T),...
                    obj(i).RootMarkerVector);
                % Sort number of time steps from largest to smallest
                b = sort(noTimeStepsVector,'descend');
                % The time steps to filter markers corresponds to the
                % noMarkers-th element of the sorted vector
                obj(i).NoTimeStepsFilter = b(noMarkers);
            end
            % Iterate through the elements of the StbStructureRun vector
            % that have less or same markers than the parameter specified
            for i=find(noMarkers>=arrayfun(@(x)...
                    length(x.RootMarkerVector),obj))'
                obj(i).NoTimeStepsFilter = length(obj(i).RootMarkerVector);
            end
        end
        
        %% Phase averaging method
        function phaseAveragedStruct = phaseAveragedDisplacements(...
                obj,noPhases,referencePeriod,timeOffset,noMarkers)
            %phaseAveragedStructure Perform phase average of markers
            %displacement
            %   phaseAveragedStructure = phaseAveragedDisplacements(...
            %   obj,noPhases,referencePeriod,timeOffset,noMarkers)
            fprintf('Phase average of structural displacements\n')
            % If time offset is not specified then assume zero time offset
            % for all STB runs
            if nargin<4
                timeOffset = zeros(size(obj));
            end
            % If number of markers is not specified
            if nargin<5
                % Use a number of time steps equal to the average value
                % plus a standard deviation among all markers to filter out
                % markers
                noTimeStepsArray = arrayfun(@(x) mean(arrayfun(@(y)...
                    length(y.T),x.RootMarkerVector))+std(arrayfun(@(y)...
                    length(y.T),x.RootMarkerVector)),obj,...
                    'UniformOutput',false);
                [obj.NoTimeStepsFilter] = noTimeStepsArray{:};
            else
                % If number of markers is specified call the
                % filterMarkersByNo method to set the number of time steps
                % to filter out markers
                obj.filterMarkersByNo(noMarkers);
            end
            % Define the nondimensional time vector corresponding to the
            % number of phases desired
            nondimensionalTimeVector = 0:1/(noPhases):1;
            % Iterate through the StbRun objects to define the
            % nondimensional time instants of STB particles
            for i=length(obj):-1:1
                fprintf('Preprocessing STB run [%d/%d]\n',length(obj)-i+1,...
                    length(obj))
                % Initialize cell arrays
                singleRunMarkerNondimensionalTimeArray =...
                    cell(size(obj(i).MarkerVector));
                singleRunMarkerXarray = cell(size(obj(i).MarkerVector));
                singleRunMarkerYarray = cell(size(obj(i).MarkerVector));
                singleRunMarkerZarray = cell(size(obj(i).MarkerVector));
                % Define the time interval for the calculation of the
                % phase-averaged displacements
                fitTimeInterval = referencePeriod(i)/4/2;
                % Iterate through the markers of current StbStructureRun
                % object
                for j=length(obj(i).MarkerVector):-1:1
                    fprintf('-> marker [%d/%d]\n',...
                        length(obj(i).MarkerVector)-j+1,...
                        length(obj(i).MarkerVector))
                    % Assemble a vector with the measurement instants of
                    % the particles of the current StbRun object,
                    % zeroing such instants with the instant of the first
                    % phase of the reference periodic signal
                    singleRunMarkerTimeArray = obj(i).MarkerVector(j).T+...
                        timeOffset(i);
                    % Nondimensionalize the time instants with the cycle
                    % period and store result in a cell array where each
                    % cell corresponds to a different StbRun object
                    singleRunMarkerNondimensionalTimeArray{j,1} =...
                        singleRunMarkerTimeArray-floor(...
                        singleRunMarkerTimeArray/referencePeriod(i))*...
                        referencePeriod(i);
                    % Find indices of particles located within half
                    % temporal bin from the start of the average cycle.
                    % These particles have to be repeated at the end of the
                    % average cycle so that temporal bins centred towards
                    % the end of the cycle can include such particles
                    appendLogicalVector =...
                        singleRunMarkerNondimensionalTimeArray{j,1}<=...
                        fitTimeInterval;
                    % Find indices of particles located within half
                    % temporal bin from the end of the average cycle. These
                    % particles have to be repeated at the start of the
                    % average cycle so that temporal bins centred towards
                    % the start of the cycle can include such particles
                    prependLogicalVector = 1-...
                        singleRunMarkerNondimensionalTimeArray{j,1}<=...
                        fitTimeInterval;
                    % Repeat positional information of identified markers
                    singleRunMarkerXarray{j,1} = [...
                        obj(i).MarkerVector(j).X(...
                        prependLogicalVector);obj(i).MarkerVector(j).X;...
                        obj(i).MarkerVector(j).X(appendLogicalVector)];
                    singleRunMarkerYarray{j,1} = [...
                        obj(i).MarkerVector(j).Y(...
                        prependLogicalVector);obj(i).MarkerVector(j).Y;...
                        obj(i).MarkerVector(j).Y(appendLogicalVector)];
                    singleRunMarkerZarray{j,1} = [...
                        obj(i).MarkerVector(j).Z(...
                        prependLogicalVector);obj(i).MarkerVector(j).Z;...
                        obj(i).MarkerVector(j).Z(appendLogicalVector)];
                    % Repeat temporal information of identified markers
                    singleRunMarkerNondimensionalTimeArray{j,1} = [...
                        singleRunMarkerNondimensionalTimeArray{j,1}(...
                        prependLogicalVector)-1;...
                        singleRunMarkerNondimensionalTimeArray{j,1};...
                        singleRunMarkerNondimensionalTimeArray{j,1}(...
                        appendLogicalVector)+1];
                end
                % Assemble cell arrays with a number of cells equal to the
                % number of StbStructureRun objects. Each cell contains a
                % cella array with the temporal and the positional
                % information of the markers
                completeMarkerNondimensionalTimeArray{i,1} =...
                    singleRunMarkerNondimensionalTimeArray;
                completeMarkerXarray{i,1} = singleRunMarkerXarray;
                completeMarkerYarray{i,1} = singleRunMarkerYarray;
                completeMarkerZarray{i,1} = singleRunMarkerZarray;
            end
            % Generate cell arrays with the temporal and positional
            % information of all markers, concatenating the cell arrays
            % corresponding to different StbStructureRun objects. The cell
            % arrays with information of all markers are ordered in
            % descending order in terms of number of time steps
            completeMarkerNondimensionalTimeArray = vertcat(...
                completeMarkerNondimensionalTimeArray{:});
            [~,I] = sort(cellfun(@(x) length(x),...
                completeMarkerNondimensionalTimeArray),'descend');
            completeMarkerNondimensionalTimeArray =...
                completeMarkerNondimensionalTimeArray(I);
            completeMarkerXarray = vertcat(completeMarkerXarray{:});
            completeMarkerXarray = completeMarkerXarray(I);
            completeMarkerYarray = vertcat(completeMarkerYarray{:});
            completeMarkerYarray = completeMarkerYarray(I);
            completeMarkerZarray = vertcat(completeMarkerZarray{:});
            completeMarkerZarray = completeMarkerZarray(I);
            %
%             % Plot markers before phase averaging (debugging)
%             figure
%             hold on
%             arrayfun(@(x) scatter3(completeMarkerXarray{x},...
%                 completeMarkerYarray{x},...
%                 completeMarkerZarray{x},'MarkerEdgeColor',...
%                 'k','MarkerFaceColor','b'),...
%                 1:length(completeMarkerZarray));
%             axis image
            %
            % Carry out displacement phase average by means of polynomial
            % fitting
            [markerAveragedXarray,markerAveragedYarray,...
                markerAveragedZarray] = polynomialPhaseAverage(...
                completeMarkerNondimensionalTimeArray,noPhases,...
                nondimensionalTimeVector,fitTimeInterval,...
                completeMarkerXarray,completeMarkerYarray,...
                completeMarkerZarray);
            %
%             % Plot markers after phase averaging (debugging)
%             figure
%             hold on
%             arrayfun(@(x) scatter3(markerAveragedXarray{x},...
%                 markerAveragedYarray{x},...
%                 markerAveragedZarray{x},'MarkerEdgeColor',...
%                 'k','MarkerFaceColor','b'),...
%                 1:length(markerAveragedXarray));
%             axis image
            %
            % Set starting index for analysis of coincident markers
            i = 1;
            % Initialize the vector of residual indices
            residualIndexVector = (1:length(markerAveragedXarray))';
            % Initialize the vector of indices to be eliminated from
            % the cell arrays used for the phase-average
            eliminateIndexVector = [];
            % Iterate through the markers to find coincident ones
            fprintf('Merging markers and averaging positions\n')
            while ~isempty(residualIndexVector)
                % If current marker is whithin the markers still to be
                % checked
                if any(residualIndexVector==i)
                    % Eliminate index of analyzed marker from the vector of
                    % residual indices
                    residualIndexVector(residualIndexVector==i) = [];
                    % Calculate a cell array containing the distance of
                    % the current marker from the remaining markers at
                    % all phases
                    distanceArray = arrayfun(@(x) vecnorm(...
                        [markerAveragedXarray{x},...
                        markerAveragedYarray{x},...
                        markerAveragedZarray{x}]-...
                        [markerAveragedXarray{i},...
                        markerAveragedYarray{i},...
                        markerAveragedZarray{i}],2,2),...
                        residualIndexVector,'UniformOutput',false);
                    % Find markers with a minimum distance from the
                    % current marker equal or lower than the marker
                    % search radius
                    sameMarkerLogical = cellfun(@(x) mean(x)+std(x),...
                        distanceArray)<=obj(1).SearchRadius;
                    if any(sameMarkerLogical)
                        % If at least one marker to be merged with current
                        % analyzed marker is found
                        %
%                         % Plot analyzed marker (debugging)
%                         scatter3(markerAveragedXarray{i},...
%                             markerAveragedYarray{i},...
%                             markerAveragedZarray{i},...
%                             'MarkerEdgeColor','k','MarkerFaceColor','g')
%                         drawnow
%                         % Plot found markers to be merged
%                         arrayfun(@(x) scatter3(markerAveragedXarray{x},...
%                             markerAveragedYarray{x},...
%                             markerAveragedZarray{x},'MarkerEdgeColor',...
%                             'k','MarkerFaceColor','r'),...
%                             residualIndexVector(sameMarkerLogical));
%                         drawnow
                        %
                        % Combine temporal and positional data of the found
                        % markers with the one of the current analyzed
                        % marker
                        completeMarkerNondimensionalTimeArray{i} = [...
                            completeMarkerNondimensionalTimeArray{i};...
                            vertcat(...
                            completeMarkerNondimensionalTimeArray{...
                            residualIndexVector(sameMarkerLogical)})];
                        completeMarkerXarray{i} = [...
                            completeMarkerXarray{i};vertcat(...
                            completeMarkerXarray{residualIndexVector(...
                            sameMarkerLogical)})];
                        completeMarkerYarray{i} = [...
                            completeMarkerYarray{i};vertcat(...
                            completeMarkerYarray{residualIndexVector(...
                            sameMarkerLogical)})];
                        completeMarkerZarray{i} = [...
                            completeMarkerZarray{i};vertcat(...
                            completeMarkerZarray{residualIndexVector(...
                            sameMarkerLogical)})];
                        % Update vector for the elimination of duplicate
                        % markers
                        eliminateIndexVector = [eliminateIndexVector;...
                            residualIndexVector(sameMarkerLogical)];
                        % Update vector with the indices of the markers
                        % that have yet to be checked
                        residualIndexVector = residualIndexVector(...
                            ~sameMarkerLogical);
                        fprintf(['-> merged %d markers, remaining %d ',...
                            'markers to check\n'],length(find(...
                            sameMarkerLogical))+1,length(...
                            residualIndexVector))
                    end
                end
                % Increase marker index to analyze next marker
                i = i+1;
            end
            % Eliminate cells corresponding to merged markers
            completeMarkerNondimensionalTimeArray(...
                eliminateIndexVector) = [];
            completeMarkerXarray(eliminateIndexVector) = [];
            completeMarkerYarray(eliminateIndexVector) = [];
            completeMarkerZarray(eliminateIndexVector) = [];
            % Perform final phase average on the reduced set of markers
            [markerAveragedXarray,markerAveragedYarray,...
                markerAveragedZarray,markerAveragedXvelocityArray,...
                markerAveragedYvelocityArray,...
                markerAveragedZvelocityArray,...
                markerAveragedXaccelerationArray,...
                markerAveragedYaccelerationArray,...
                markerAveragedZaccelerationArray] =...
                polynomialPhaseAverage(...
                completeMarkerNondimensionalTimeArray,noPhases,...
                nondimensionalTimeVector,fitTimeInterval,...
                completeMarkerXarray,completeMarkerYarray,...
                completeMarkerZarray);
            % Filter out markers where NaN is present in at least one
            % phase
            fprintf('Filter out markers with NaN in at least one phase\n')
            filterVector = cellfun(@(x) ~any(isnan(x)),...
                markerAveragedXarray);
            markerAveragedXarray = markerAveragedXarray(filterVector);
            markerAveragedYarray = markerAveragedYarray(filterVector);
            markerAveragedZarray = markerAveragedZarray(filterVector);
            markerAveragedXvelocityArray =...
                markerAveragedXvelocityArray(filterVector);
            markerAveragedYvelocityArray =...
                markerAveragedYvelocityArray(filterVector);
            markerAveragedZvelocityArray =...
                markerAveragedZvelocityArray(filterVector);
            markerAveragedXaccelerationArray =...
                markerAveragedXaccelerationArray(filterVector);
            markerAveragedYaccelerationArray =...
                markerAveragedYaccelerationArray(filterVector);
            markerAveragedZaccelerationArray =...
                markerAveragedZaccelerationArray(filterVector);
            %
%             % Plot markers after final phase average (debugging)
%             figure
%             hold on
%             arrayfun(@(x) scatter3(markerAveragedXarray{x},...
%                 markerAveragedYarray{x},markerAveragedZarray{x},...
%                 'MarkerEdgeColor','k','MarkerFaceColor','b'),...
%                 1:length(markerAveragedXarray));
%             axis image
            %
            % Rearrange the data in a cell array with a number of rows
            % equal to the number of phases and having in each element
            % a vector with the data of all markers
            fprintf('Assemble output structure\n')
            for i=noPhases+1:-1:1
                for j=length(markerAveragedXarray):-1:1
                    markerAveragedXoutArray{i,1}(j,1) =...
                        markerAveragedXarray{j}(i);
                    markerAveragedYoutArray{i,1}(j,1) =...
                        markerAveragedYarray{j}(i);
                    markerAveragedZoutArray{i,1}(j,1) =...
                        markerAveragedZarray{j}(i);
                    markerAveragedXvelocityOutArray{i,1}(j,1) =...
                        markerAveragedXvelocityArray{j}(i);
                    markerAveragedYvelocityOutArray{i,1}(j,1) =...
                        markerAveragedYvelocityArray{j}(i);
                    markerAveragedZvelocityOutArray{i,1}(j,1) =...
                        markerAveragedZvelocityArray{j}(i);
                    markerAveragedXaccelerationOutArray{i,1}(j,1) =...
                        markerAveragedXaccelerationArray{j}(i);
                    markerAveragedYaccelerationOutArray{i,1}(j,1) =...
                        markerAveragedYaccelerationArray{j}(i);
                    markerAveragedZaccelerationOutArray{i,1}(j,1) =...
                        markerAveragedZaccelerationArray{j}(i);
                end
            end
            % Assemble output structure
            phaseAveragedStruct = struct(...
                'xPosition',markerAveragedXoutArray,...
                'yPosition',markerAveragedYoutArray,...
                'zPosition',markerAveragedZoutArray,...
                'xVelocity',markerAveragedXvelocityOutArray,...
                'yVelocity',markerAveragedYvelocityOutArray,...
                'zVelocity',markerAveragedZvelocityOutArray,...
                'xAcceleration',markerAveragedXaccelerationOutArray,...
                'yAcceleration',markerAveragedYaccelerationOutArray,...
                'zAcceleration',markerAveragedZaccelerationOutArray);
        end
        
        %% Time averaging method
        function timeAveragedStruct = timeAveragedDisplacements(...
                obj,noMarkers)
            %phaseAveragedStructure Perform phase average of markers
            %displacement
            %   phaseAveragedStructure = phaseAveragedDisplacements(...
            %   obj,noPhases,referencePeriod,timeOffset,noMarkers)
            fprintf('Time average of structural displacements\n')
            % If number of markers is not specified
            if nargin>=2
                % If number of markers is specified call the
                % filterMarkersByNo method to set the number of time steps
                % to filter out markers
                obj.filterMarkersByNo(noMarkers);
            end
            % Iterate through the StbStructureRun objects
            for i=length(obj):-1:1
                fprintf('Preprocessing STB run [%d/%d]\n',...
                    length(obj)-i+1,length(obj))
                % Initialize cell arrays
                singleRunMarkerXarray = cell(size(obj(i).MarkerVector));
                singleRunMarkerYarray = cell(size(obj(i).MarkerVector));
                singleRunMarkerZarray = cell(size(obj(i).MarkerVector));
                % Iterate through the markers of current StbStructureRun
                % object
                for j=length(obj(i).MarkerVector):-1:1
                    fprintf('-> marker [%d/%d]\n',...
                        length(obj(i).MarkerVector)-j+1,...
                        length(obj(i).MarkerVector))
                    % Assemble cell arrays with the positional information
                    % of each marker
                    singleRunMarkerXarray{j,1} = obj(i).MarkerVector(j).X;
                    singleRunMarkerYarray{j,1} = obj(i).MarkerVector(j).Y;
                    singleRunMarkerZarray{j,1} = obj(i).MarkerVector(j).Z;
                end
                % Assemble cell arrays with positional information of all
                % markers from all StbStructureRun objects
                completeMarkerXarray{i,1} = singleRunMarkerXarray;
                completeMarkerYarray{i,1} = singleRunMarkerYarray;
                completeMarkerZarray{i,1} = singleRunMarkerZarray;
            end
            % Generate cell arrays with positional information of all
            % markers, concatenating the cell arrays corresponding to
            % different StbStructureRun objects. The cell arrays with
            % information of all markers are ordered in descending order in
            % terms of number of time steps
            completeMarkerXarray = vertcat(completeMarkerXarray{:});
            [~,I] = sort(cellfun(@(x) length(x),completeMarkerXarray),...
                'descend');
            completeMarkerXarray = completeMarkerXarray(I);
            completeMarkerYarray = vertcat(completeMarkerYarray{:});
            completeMarkerYarray = completeMarkerYarray(I);
            completeMarkerZarray = vertcat(completeMarkerZarray{:});
            completeMarkerZarray = completeMarkerZarray(I);
            % Carry out displacement time average
            for i=length(completeMarkerXarray):-1:1
                markerAveragedXvector(i) = mean(completeMarkerXarray{i});
                markerAveragedYvector(i) = mean(completeMarkerYarray{i});
                markerAveragedZvector(i) = mean(completeMarkerZarray{i});
            end
            % Set starting index for analysis of coincident markers
            i = 1;
            % Initialize the vector of residual indices
            residualIndexVector = (1:length(markerAveragedXvector))';
            % Initialize the vector of indices to be eliminated from
            % the cell arrays used for the phase-average
            eliminateIndexVector = [];
            % Iterate through the markers to find coincident ones
            fprintf('Merging markers and averaging positions\n')
            while ~isempty(residualIndexVector)
                % If current marker is whithin the markers still to be
                % checked
                if any(residualIndexVector==i)
                    % Eliminate index of analyzed marker from the vector of
                    % residual indices
                    residualIndexVector(residualIndexVector==i) = [];
                    % Calculate a cell array containing the distance of
                    % the current marker from the remaining markers at
                    % all phases
                    distanceArray = arrayfun(@(x) vecnorm(...
                        [markerAveragedXvector(x),...
                        markerAveragedYvector(x),...
                        markerAveragedZvector(x)]-...
                        [markerAveragedXvector(i),...
                        markerAveragedYvector(i),...
                        markerAveragedZvector(i)],2,2),...
                        residualIndexVector,'UniformOutput',false);
                    % Find markers with a minimum distance from the
                    % current marker equal or lower than the marker
                    % search radius
                    sameMarkerLogical = cellfun(@(x) mean(x)+std(x),...
                        distanceArray)<=obj(1).SearchRadius;
                    if any(sameMarkerLogical)
                        % If at least one marker to be merged with current
                        % analyzed marker is found, combine the positional
                        % data of the found markers with the one of the 
                        % current analyzed marker
                        completeMarkerXarray{i} = [...
                            completeMarkerXarray{i};vertcat(...
                            completeMarkerXarray{residualIndexVector(...
                            sameMarkerLogical)})];
                        completeMarkerYarray{i} = [...
                            completeMarkerYarray{i};vertcat(...
                            completeMarkerYarray{residualIndexVector(...
                            sameMarkerLogical)})];
                        completeMarkerZarray{i} = [...
                            completeMarkerZarray{i};vertcat(...
                            completeMarkerZarray{residualIndexVector(...
                            sameMarkerLogical)})];
                        % Update vector for the elimination of duplicate
                        % markers
                        eliminateIndexVector = [eliminateIndexVector;...
                            residualIndexVector(sameMarkerLogical)];
                        % Update vector with the indices of the markers
                        % that have yet to be checked
                        residualIndexVector = residualIndexVector(...
                            ~sameMarkerLogical);
                        fprintf(['-> merged %d markers, remaining %d ',...
                            'markers to check\n'],length(find(...
                            sameMarkerLogical))+1,length(...
                            residualIndexVector))
                    end
                end
                % Increase marker index to analyze next marker
                i = i+1;
            end
            % Eliminate cells corresponding to merged markers
            completeMarkerXarray(eliminateIndexVector) = [];
            completeMarkerYarray(eliminateIndexVector) = [];
            completeMarkerZarray(eliminateIndexVector) = [];
            % Perform final time average on the reduced set of markers
            for i=length(completeMarkerXarray):-1:1
                markerAveragedXvector(i) = mean(completeMarkerXarray{i});
                markerAveragedYvector(i) = mean(completeMarkerYarray{i});
                markerAveragedZvector(i) = mean(completeMarkerZarray{i});
            end
            % Filter out markers where NaN is present in at least one
            % phase
            fprintf('Filter out markers with NaN in at least one phase\n')
            filterVector = ~isnan(markerAveragedXvector);
            markerAveragedXvector = markerAveragedXvector(filterVector);
            markerAveragedYvector = markerAveragedYvector(filterVector);
            markerAveragedZvector = markerAveragedZvector(filterVector);
            % Assemble output structure
            timeAveragedStruct = struct(...
                'xPosition',markerAveragedXvector,...
                'yPosition',markerAveragedYvector,...
                'zPosition',markerAveragedZvector);
        end
        
        %% plotMarkers method
        function h = plotMarkers(obj,varargin)
            markerVector = vertcat(obj.MarkerVector);
            h = markerVector.plot(varargin{:});
        end
        
        %% saveobj method
        function s = saveobj(obj)
            fprintf('Saving StbRun object\n')
            % Store path to original .dat file
            fprintf('-> path to original .dat file\n')
            s.rootFile = obj.RootFile;
            % Store an array with all particles' positional and velocity
            % information
            fprintf(['-> array with particles'' positional and ',...
                'velocity information\n'])
            s.particleDataArray = [vertcat(obj.ParticleVector.X),...
                vertcat(obj.ParticleVector.Y),...
                vertcat(obj.ParticleVector.Z),...
                vertcat(obj.ParticleVector.I),...
                vertcat(obj.ParticleVector.U),...
                vertcat(obj.ParticleVector.V),...
                vertcat(obj.ParticleVector.W),...
                vertcat(obj.ParticleVector.VelocityMagnitude),...
                vertcat(obj.ParticleVector.Ax),...
                vertcat(obj.ParticleVector.Ay),...
                vertcat(obj.ParticleVector.Az),...
                vertcat(obj.ParticleVector.AccelerationMagnitude)];
            % Store a vector with the time instants of all particles
            fprintf('-> array with particles'' time instant\n')
            s.timeVector = vertcat(obj.ParticleVector.T);
            % Store a vector with the time steps of all particles
            fprintf('-> array with particles'' time step number\n')
            s.timeStepNoVector = vertcat(obj.ParticleVector.TimeStepNo);
            % Store a vector with the track id of all particles
            fprintf('-> array with particles'' track id\n')
            s.trackIdVector = vertcat(obj.ParticleVector.TrackId);
            % Store a cell array with the track id of all markers
            fprintf('-> marker info\n')
            s.rootMarkerTrackIdArray = arrayfun(@(x)...
                [x.TrackVector.TrackId],obj.RootMarkerVector,...
                'UniformOutput',false);
            % Store additional parameters of StbStructureRun object
            s.searchRadius = obj.SearchRadius;
            s.fitTimeInterval = obj.FitTimeInterval;
            s.noTimeStepsFilter = obj.NoTimeStepsFilter;
        end
    end
    
    %% loadobj method
    methods(Static)
        function obj = loadobj(s)
            newObj = StbStructureRun;
            % Regenerate Particle object
            fprintf('Regenerating array of StbParticle objects\n')
            particleVector = StbParticle(...
                s.particleDataArray(:,1),...
                s.particleDataArray(:,2),...
                s.particleDataArray(:,3),...
                s.particleDataArray(:,4),...
                s.particleDataArray(:,5),...
                s.particleDataArray(:,6),...
                s.particleDataArray(:,7),...
                s.particleDataArray(:,8),...
                s.particleDataArray(:,9),...
                s.particleDataArray(:,10),...
                s.particleDataArray(:,11),...
                s.particleDataArray(:,12));
            % Assemble Particle objects per time step
            [particle2TimeStepIndexVector,uniqueTimeStepVector] =...
                findgroups(s.timeStepNoVector);
            particlePerTimeStepArray = splitapply(@(x){particleVector(x)...
                },(1:length(s.timeStepNoVector))',...
                particle2TimeStepIndexVector);
            % Regenerate TimeStep object
            fprintf('Regenerating array of StbTimeStep objects\n')
            timeStepObjVector = StbTimeStep(uniqueTimeStepVector,...
                unique(s.timeVector),particlePerTimeStepArray);
            % Assemble Particle objects per track id
            % trackIdVector includes the track id of all particles, as a consequence
            % particles from the same track generate duplicates in the vector. Once the
            % vector is obtained, it is necessary to map the indices of same track id
            % to the correct track in order to assemble the Track object with the
            % correct children Particle objects
            [particle2TrackIndexVector,uniqueTrackIdVector] =...
                findgroups(s.trackIdVector);
            particlePerTrackArray = splitapply(@(x){particleVector(x)},...
                (1:length(s.trackIdVector))',particle2TrackIndexVector);
            % Regenerate StbTrack object
            fprintf('Regenerating array of StbTrack objects\n')
            trackVector = StbTrack(num2cell(uniqueTrackIdVector),...
                particlePerTrackArray);
            % Regenerate StructuralMarker object
            fprintf('Regenerating array of StructuralMarker objects\n')
            markerTrackArray = cellfun(@(x) trackVector(x+1),...
                s.rootMarkerTrackIdArray,'UniformOutput',false);
            markerVector(length(markerTrackArray),1) = StructuralMarker;
            [markerVector.TrackVector] = markerTrackArray{:};
            % Assign generated objects to properties
            newObj.RootFile = s.rootFile;
            newObj.ParticleVector = particleVector;
            newObj.TimeStepVector = timeStepObjVector;
            newObj.TrackVector = trackVector;
            newObj.RootMarkerVector = markerVector;
            newObj.MarkerVectorProxy = newObj.RootMarkerVector;
            newObj.SearchRadiusProxy = s.searchRadius;
            newObj.FitTimeIntervalProxy = s.fitTimeInterval;
            if ~isempty(s.noTimeStepsFilter)
                newObj.NoTimeStepsFilter = s.noTimeStepsFilter;
            end
            obj = newObj;
        end
    end
end

%% Phase average function
function [markerAveragedXarray,markerAveragedYarray,...
    markerAveragedZarray,markerAveragedXvelocityArray,...
    markerAveragedYvelocityArray,markerAveragedZvelocityArray,...
    markerAveragedXaccelerationArray,markerAveragedYaccelerationArray,...
    markerAveragedZaccelerationArray] = polynomialPhaseAverage(...
    markerNondimensionalTimeArray,noPhases,nondimensionalTimeVector,...
    fitTimeInterval,markerXarray,markerYarray,markerZarray)
%polynomialPhaseAverage Perform phase average of markers' displacement by
%means of a polynomial fit
%   [markerAveragedXarray,markerAveragedYarray,...
%   markerAveragedZarray] = polynomialPhaseAverage(...
%   markerNondimensionalTimeArray,noPhases,nondimensionalTimeVector,...
%   fitTimeInterval,markerXarray,markerYarray,markerZarray)
% Iterate through the markers
fprintf('Phase average of markers\n')
parfor i=1:length(markerNondimensionalTimeArray)
    fprintf('-> [%d/%d]\n',i,length(...
        markerNondimensionalTimeArray))
    % Initialize vector of the averaged x, y and z coordinates
    markerAveragedXvector = NaN(noPhases+1,1);
    markerAveragedYvector = NaN(noPhases+1,1);
    markerAveragedZvector = NaN(noPhases+1,1);
    markerAveragedXvelocityVector = NaN(noPhases+1,1);
    markerAveragedYvelocityVector = NaN(noPhases+1,1);
    markerAveragedZvelocityVector = NaN(noPhases,1);
    markerAveragedXaccelerationVector = NaN(noPhases+1,1);
    markerAveragedYaccelerationVector = NaN(noPhases+1,1);
    markerAveragedZaccelerationVector = NaN(noPhases+1,1);
    % Iterate through the the nondimensional time vector corresponding to
    % the number of phases
    for t=1:length(nondimensionalTimeVector)
        % Find logical vector used to select the time kernel of the
        % polynomial fit
        fitLogical =...
            markerNondimensionalTimeArray{i}>=...
            nondimensionalTimeVector(t)-fitTimeInterval...
            & markerNondimensionalTimeArray{i}<=...
            nondimensionalTimeVector(t)+fitTimeInterval;
        % If at least 6 particles are present in the time interval
        % selected and if the nondimensional time instant corresponding to
        % the current phase is within the minimum and the maximum time
        % instant of the current marker, then phase average the
        % displacements
        if length(find(fitLogical))>=6 &&...
                nondimensionalTimeVector(t)>=...
                min(markerNondimensionalTimeArray{i}...
                ) && nondimensionalTimeVector(t)<=...
                max(markerNondimensionalTimeArray{i})
            % x displacement
            px = polyfit(...
                markerNondimensionalTimeArray{i...
                }(fitLogical),markerXarray{i}(...
                fitLogical),2);
            markerAveragedXvector(t) = polyval(px,...
                nondimensionalTimeVector(t));
            % y displacement
            py = polyfit(...
                markerNondimensionalTimeArray{i...
                }(fitLogical),markerYarray{i}(...
                fitLogical),2);
            markerAveragedYvector(t) = polyval(py,...
                nondimensionalTimeVector(t));
            % z displacement
            pz = polyfit(...
                markerNondimensionalTimeArray{i...
                }(fitLogical),markerZarray{i}(...
                fitLogical),2);
            markerAveragedZvector(t) = polyval(pz,...
                nondimensionalTimeVector(t));
            % x velocity
            markerAveragedXvelocityVector(t) = 2*px(1)*...
                nondimensionalTimeVector(t)+px(2);
            % x acceleration
            markerAveragedXaccelerationVector(t) = 2*px(1);
            % y velocity
            markerAveragedYvelocityVector(t) = 2*py(1)*...
                nondimensionalTimeVector(t)+py(2);
            % y acceleration
            markerAveragedYaccelerationVector(t) = 2*py(1);
            % z velocity
            markerAveragedZvelocityVector(t) = 2*pz(1)*...
                nondimensionalTimeVector(t)+pz(2);
            % z acceleration
            markerAveragedZaccelerationVector(t) = 2*pz(1);
        end
    end
    % Assign averaged positional vectors to output cell arrays
    markerAveragedXarray{i,1} = markerAveragedXvector;
    markerAveragedYarray{i,1} = markerAveragedYvector;
    markerAveragedZarray{i,1} = markerAveragedZvector;
    % Assign averaged velocities and accelerations vector to output cell
    % arrays
    markerAveragedXvelocityArray{i,1} = markerAveragedXvelocityVector;
    markerAveragedYvelocityArray{i,1} = markerAveragedYvelocityVector;
    markerAveragedZvelocityArray{i,1} = markerAveragedZvelocityVector;
    markerAveragedXaccelerationArray{i,1} =...
        markerAveragedXaccelerationVector;
    markerAveragedYaccelerationArray{i,1} =...
        markerAveragedYaccelerationVector;
    markerAveragedZaccelerationArray{i,1} =...
        markerAveragedZaccelerationVector;
end
end

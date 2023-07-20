classdef StructuralMarker < matlab.mixin.Copyable
    %Track Summary of this class goes here
    %   Detailed explanation goes here
    
    %% Properties
    properties (Dependent)
        T
        TimeStepNo
        X
        Y
        Z
        U
        V
        W
        TrackVector
    end
    properties (Hidden,SetAccess=private)
        Tproxy
        TimeStepNoProxy
        Xproxy
        Yproxy
        Zproxy
        Uproxy
        Vproxy
        Wproxy
        TrackVectorProxy
    end
    
    methods
        %% Constructor
        function obj = StructuralMarker(timeStepArray,...
                markerFitTimeInterval,searchRadius)
            %StructuralMarker Construct an instance of this class
            %   obj = StructuralMarker(timeStepArray,...
            %   markerFitTimeInterval,searchRadius)
            % If number of input arguments is not zero then initialize the
            % object array with the size of the first input
            if nargin ~= 0
                % Get tracks starting in first time step
                firstTimeStepTrackArray =...
                    [timeStepArray(1).ParticleVector.Track];
                fistTimeStepTrackCellArray =...
                    num2cell(firstTimeStepTrackArray);
                % Get tracks in the entire TimeStep vector
                residualTrackArray = timeStepArray.getTrackVector;
                % Initialize object using the number of tracks in first
                % time step
                L = length(firstTimeStepTrackArray);
                fprintf(['Generating %d markers from the tracks found ',...
                    'in the first time step\n'],L)
                obj(L,1) = StructuralMarker;
                % Assign the tracks of the first time step to the markers
                [obj.TrackVector] = fistTimeStepTrackCellArray{:};
                % Update residual tracks
                residualTrackArray = residualTrackArray(~ismembc(...
                    [residualTrackArray.TrackId],...
                    [firstTimeStepTrackArray.TrackId]));
                %
                %                 % Plot initial markers (debugging)
                %                 obj.plot
                %                 hold on
                %
                % Set badly conditioned polynomial warning (polyfit) to
                % temporarily issue errors (exceptions)
                s = warning('error',...
                    'MATLAB:polyfit:RepeatedPointsOrRescale');
                % Iterate through the time steps
                for i=2:length(timeStepArray)
                    % Print info on current time step and residual tracks
                    fprintf('Time step %d/%d\n# residual tracks: %d\n',...
                        i,length(timeStepArray),...
                        length(residualTrackArray));
                    % Retrieve tracks from current time step
                    currentTimeStepTrackArray =...
                        vertcat(timeStepArray(i).ParticleVector.Track);
                    % Find residual tracks (unassigned to any marker) in
                    % the current time steps
                    newTrackArray = residualTrackArray(ismembc(...
                        [residualTrackArray.TrackId],...
                        [currentTimeStepTrackArray.TrackId]));
                    % If new tracks are found
                    if ~isempty(newTrackArray)
                        % Find last time step of each marker
                        for j=length(obj):-1:1
                            markerLastTimeStepVector(j) = obj(j).T(end);
                        end
                        % Select only the markers with last time step
                        % smaller than the current time step and that
                        % result within the maximum track separation time
                        % allowed
                        markerSelectionLogical = timeStepArray(i).T>...
                            markerLastTimeStepVector &...
                            timeStepArray(i).T-markerLastTimeStepVector...
                            <=markerFitTimeInterval;
                        markerSelection = obj(markerSelectionLogical);
                        % If at least one marker is found, try to assign
                        % one of the unassigned tracks to one of the
                        % markers
                        if ~isempty(markerSelection)
                            [newTrackArray,residualTrackArray] =...
                                findNewTrack4Marker(markerSelection,...
                                newTrackArray,residualTrackArray,...
                                markerFitTimeInterval,searchRadius);
                            %
                            %                         else
                            %                             % Plot new markers (debugging)
                            %                             newTrackArray.plot('color','r')
                            %                             drawnow
                            %
                        end
                    end
                    if ~isempty(newTrackArray)
                        % If at the end of the iteration through the
                        % markers the vector containing the new tracks
                        % is not empty, then add a new marker
                        for j=1:length(newTrackArray)
                            fprintf(['New marker found. Current # ',...
                                'markers: %d\n'],length(obj)+1)
                            obj(end+1).TrackVector = newTrackArray(j);
                            % Update residualTrackArray
                            residualTrackArray(residualTrackArray==...
                                newTrackArray(j)) = [];
                        end
                    end
                end
                % Restore the warning back to their previous (non-error)
                % state
                warning(s);
            end
        end
        
        %% Filter markers by tracks number
        function filteredMarkerVector = filterMarkersByNoTimeSteps(obj,...
                noTimeSteps)
            fprintf('Finding markers with at least %d time steps\n',...
                noTimeSteps)
            filteredMarkerVector = obj(arrayfun(@(x)...
                length(x.TimeStepNo)>=noTimeSteps,obj));
        end
        
        %% TrackVector set and get method
        function set.TrackVector(obj,trackVector)
            % Assign TrackVectorProxy property
            obj.TrackVectorProxy = trackVector;
            % Assign proxy properties
            obj.Tproxy = vertcat(trackVector.T);
            obj.TimeStepNoProxy = vertcat(trackVector.TimeStepNo);
            obj.Xproxy = vertcat(trackVector.X);
            obj.Yproxy = vertcat(trackVector.Y);
            obj.Zproxy = vertcat(trackVector.Z);
            obj.Uproxy = vertcat(trackVector.U);
            obj.Vproxy = vertcat(trackVector.V);
            obj.Wproxy = vertcat(trackVector.W);
        end
        function trackVector = get.TrackVector(obj)
            trackVector = obj.TrackVectorProxy;
        end
        
        %% Position and velocities get methods
        function t = get.T(obj)
            t = obj.Tproxy;
        end
        function timeStepNo = get.TimeStepNo(obj)
            timeStepNo = obj.TimeStepNoProxy;
        end
        function x = get.X(obj)
            x = obj.Xproxy;
        end
        function y = get.Y(obj)
            y = obj.Yproxy;
        end
        function z = get.Z(obj)
            z = obj.Zproxy;
        end
        function u = get.U(obj)
            u = obj.Uproxy;
        end
        function v = get.V(obj)
            v = obj.Vproxy;
        end
        function w = get.W(obj)
            w = obj.Wproxy;
        end
        
        %% plot method
        function h = plot(obj,varargin)
            % Create an InputParser object
            p = inputParser;
            % Add inputs to the parsing scheme
            defaultColor = [0,.75,.75];
            addRequired(p,'obj',@(obj)isa(obj,'StructuralMarker'));
            addParameter(p,'noTimeSteps2Plot',[]);
            addParameter(p,'color',defaultColor,@(x)...
                isnumeric(x)||ischar(x))
            addParameter(p,'targetAxes',gca)
            % Set properties to adjust parsing
            p.KeepUnmatched = true;
            % Parse the inputs
            parse(p,obj,varargin{:})
            % If no parameter noTimeSteps2Plot is given, plot all timesteps
            % of the markers
            if isempty(p.Results.noTimeSteps2Plot)
                timeStepIndexVector = arrayfun(@(x) 1:length(x.T),...
                    obj,'UniformOutput',false);
            else
                timeStepIndexVector = arrayfun(@(x) 1:...
                    p.Results.noTimeSteps2Plot,obj,'UniformOutput',false);
            end
            % Find coordinates corresponding to time steps to plot
            xVector = cell2mat(arrayfun(@(x) obj(x).X(...
                timeStepIndexVector{x}),(1:length(obj))',...
                'UniformOutput',false));
            yVector = cell2mat(arrayfun(@(x) obj(x).Y(...
                timeStepIndexVector{x}),(1:length(obj))',...
                'UniformOutput',false));
            zVector = cell2mat(arrayfun(@(x) obj(x).Z(...
                timeStepIndexVector{x}),(1:length(obj))',...
                'UniformOutput',false));
            % Generate scatter plot
            h = scatter3(p.Results.targetAxes,xVector,yVector,zVector,...
                'MarkerEdgeColor','k','MarkerFaceColor',p.Results.color);
            % Adjust plot
            makePlotNicer(struct('txtXlabel','$x$ [m]',...
                'txtYlabel','$y$ [m]','txtZlabel','$z$ [m]'))
            axis image
        end
        
        %% selectMarker method
        function marker = selectMarker(obj)
            h = figure;
            obj.plot('noTimeSteps2Plot',1);
            datacursormode on
            dcm_obj = datacursormode(h);
            cursorInfo = [];
            while isempty(cursorInfo)
                waitforbuttonpress;
                cursorInfo = getCursorInfo(dcm_obj);
            end
            firstParticleArray = arrayfun(@(x)...
                x.TrackVector(1).ParticleVector(1),obj);
            selectedFirstParticle = firstParticleArray(...
                [firstParticleArray.X]==cursorInfo.Position(1)&...
                [firstParticleArray.Y]==cursorInfo.Position(2)&...
                [firstParticleArray.Z]==cursorInfo.Position(3));
            marker = obj(selectedFirstParticle.Track==arrayfun(@(x)...
                x.TrackVector(1),obj));
            close(h)
        end
        
        %% Join markers method
        function marker = joinMarkers(obj,noMarkers)
            if nargin>1 && noMarkers>1
                proxyObj = copy(obj);
                for i=noMarkers:-1:1
                    fprintf('Select marker %d/%d\n',i-noMarkers+1,...
                        noMarkers)
                    marker2JoinArray(i) = proxyObj.selectMarker;
                    proxyObj(proxyObj==marker2JoinArray(i)) = [];
                end
                [~,I] = sort(vertcat(marker2JoinArray.T));
                trackArray = vertcat(marker2JoinArray.TrackVector);
                unsortedParticleArray = vertcat(trackArray.ParticleVector);
                trackArray(1).ParticleVector = unsortedParticleArray(I);
                marker2JoinArray(1).TrackVector = trackArray(1);
                marker = marker2JoinArray(1);
            end
        end
        
        %% Return time history of absolute displacement from initial point
        function [timeVector,displacementArray] = timeDisplacement(obj)
            timeVector = obj.T;
            particleXyzArray = [obj.X,obj.Y,obj.Z];
            initialPointXyz = particleXyzArray(1,:);
            displacementArray = particleXyzArray-repmat(initialPointXyz,...
                size(particleXyzArray,1),1);
        end
    end
end

function [trackVectorOut,residualTrackArray] = findNewTrack4Marker(...
    markerArray,newTrackArray,residualTrackArray,fitTimeInterval,...
    searchRadius)
% Initialize output track vector
trackVectorOut = {};
% Iterate through the markers
for i=length(markerArray):-1:1
    % Find logical vector indicating time steps included in the time
    % interval used for polynomial fit of the marker trajectory
    fitLogicalVector = markerArray(i).T>=markerArray(i).T(end)-...
        fitTimeInterval;
    % Retrieve time vector for the polynomial fit
    fitTimeVector = markerArray(i).T(fitLogicalVector);
    % Perform a second order polynomial fit in x, y and z (quadratic
    % trajectory in time). The rows of p have the coefficients of the
    % regression of each variable. If polynomial is badly conditioned, then
    % perform a first oredr polynomial fit
    try
        p{1} = polyfit(fitTimeVector,markerArray(i).X(fitLogicalVector),2);
    catch
        p{1} = polyfit(fitTimeVector,markerArray(i).X(fitLogicalVector),1);
    end
    try
        p{2} = polyfit(fitTimeVector,markerArray(i).Y(fitLogicalVector),2);
    catch
        p{2} = polyfit(fitTimeVector,markerArray(i).Y(fitLogicalVector),1);
    end
    try
        p{3} = polyfit(fitTimeVector,markerArray(i).Z(fitLogicalVector),2);
    catch
        p{3} = polyfit(fitTimeVector,markerArray(i).Z(fitLogicalVector),1);
    end
    % Iterate through the unassigned tracks
    for j=length(newTrackArray):-1:1
        % Initialize fittedXyzNewTrackTimeStepArray for current iteration
        fittedXyzNewTrackTimeStepArray = zeros(3,...
            length(newTrackArray(j).T));
        % Evaluate x, y, z of the fitted trajectory at the time steps
        % of the new track
        for k=3:-1:1
            % The columns of fittedXyzNewTrackTimeStepArray have the fitted
            % x, y and z and each row represent a different time step
            fittedXyzNewTrackTimeStepArray(k,:) = polyval(p{k},...
                newTrackArray(j).T)';
        end
        % Evaluate the mean euclidean distance between the points of
        % the new track and the ones predicted by the polynomial fit.
        % The rows of meanDistanceArray represents the unassigned tracks,
        % while the columns represent the different markers
        meanDistanceArray(j,i) = mean(vecnorm([newTrackArray(j).X,...
            newTrackArray(j).Y,newTrackArray(j).Z]'-...
            fittedXyzNewTrackTimeStepArray));
        % Evaluate standard deviation
        stdDistanceArray(j,i) = std(vecnorm([newTrackArray(j).X,...
            newTrackArray(j).Y,newTrackArray(j).Z]'-...
            fittedXyzNewTrackTimeStepArray));
    end
end
% Find minimum mean distance of each track from the predicted points of the
% markers
[minDistanceVector,markerIndexVector] = min(meanDistanceArray,[],2);
% Increase mean minimum distances with one standard deviation
for i=1:length(minDistanceVector)
    minDistanceVector(i) = minDistanceVector(i)+...
        stdDistanceArray(i,markerIndexVector(i));
end
% Find the indices of the markers having a mean distance from the
% corresponding track smaller than the search radius
validMarkerLogical = minDistanceVector<=searchRadius;
validMarkerIndexVector = markerIndexVector(validMarkerLogical);
% Find duplicate values in markerIndexVector, corresponding to multiple
% tracks having minimum mean distance from the same marker
[~,ia] = unique(validMarkerIndexVector);
% Duplicate indices
duplicateIndexVector = setdiff(1:length(validMarkerIndexVector),ia);
% Repeated markers
duplicateMarkerIndexVector = unique(validMarkerIndexVector(...
    duplicateIndexVector));
% Iterate through repeated markers
for i=length(duplicateMarkerIndexVector):-1:1
    % Find the index of the marker with the smallest distance within the
    % subset of the repeated markers
    [~,I] = min(minDistanceVector(markerIndexVector==...
        duplicateMarkerIndexVector(i)));
    % Map index of subset to index of entire marker set
    k = find(markerIndexVector==duplicateMarkerIndexVector(i));
    I = k(I);
    %
    %     % Plot track assigned to marker (debugging)
    %     newTrackArray(I).plot('color','g')
    %     drawnow
    %     % Plot track not assigned to marker
    %     newTrackArray(k(k~=I)).plot('color','r')
    %     drawnow
    %
    % Assign track to marker
    markerArray(markerIndexVector(I)).TrackVector =...
        [markerArray(markerIndexVector(I)).TrackVector;newTrackArray(I)];
    % Update array of resiudal tracks
    residualTrackArray(residualTrackArray==newTrackArray(I)) = [];
    % Update array of unassigned tracks
    trackVectorOut{i} = newTrackArray(k(k~=I));
end
% Iterate through the remaining markers having minimum mean distance from
% corresponding track that is smaller than the search radius
for i=find(minDistanceVector<=searchRadius &...
        ~ismember(markerIndexVector,duplicateMarkerIndexVector))'
    %
    %     % Plot track assigned to marker (debugging)
    %     newTrackArray(i).plot('color','g')
    %     drawnow
    %
    % Assign track to marker
    markerArray(markerIndexVector(i)).TrackVector =...
        [markerArray(markerIndexVector(i)).TrackVector;newTrackArray(i)];
    % Update array of resiudal tracks
    residualTrackArray(residualTrackArray==newTrackArray(i)) = [];
end
%
% Plot tracks with minimum mean distance from markers predicted points
% larger than search radius (debugging)
% for i=find(minDistanceVector>searchRadius)'
%     newTrackArray(i).plot('color','r')
%     drawnow
% end
%
% Update array of unassigned tracks
trackVectorOut = [vertcat(trackVectorOut{:});newTrackArray(...
    minDistanceVector>searchRadius)];
end

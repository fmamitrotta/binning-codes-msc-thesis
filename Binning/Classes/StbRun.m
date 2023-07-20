classdef StbRun < matlab.mixin.Copyable
    %StbRun Class for the handling of Shake-the-Box track data
    %   Detailed explanation goes here
    
    properties
        RootFile = '';
        ParticleVector  % vector with Particle objects
        TimeStepVector  % vector with TimeStep objects
        TrackVector     % vector with Track objects
    end
    
    methods
        %% Constructor
        function obj = StbRun(stbFilepathArray)
            %StbRun Construct an instance of this class
            %   obj = StbRun(stbFilepathArray) generates an array of StbRun
            %   objects from a cell array containing the paths to the STB
            %   .dat files
            % If number of input arguments is not zero initialize the
            % object array
            if nargin ~= 0 && ~isempty(stbFilepathArray)
                switch class(stbFilepathArray)
                    case 'cell'
                        % If stbFilepathArray is a cell array initialize
                        % object with same dimension as stbFilepathArray
                        [m,n] = size(stbFilepathArray);
                        tempCellArray = arrayfun(@(x) arrayfun(@(y)...
                            copy(obj),1:m)',1:n,'UniformOutput',false);
                        obj = horzcat(tempCellArray{:});
                    case 'char'
                        % If stbFilepathArray is a char initialize object
                        % as 1x1 array
                        m = 1;
                        n = 1;
                        stbFilepathArray = textscan(...
                            stbFilepathArray,'%s','Delimiter','\n');
                        stbFilepathArray = stbFilepathArray{:};
                    otherwise
                        % Throw error if stbFilepathArray is neither a cell
                        % nor a char
                        error(['Error. \nInput must be a cell or a ',...
                            'char, not a %s.'],class(stbFilepathArray))
                end
                % Iterate through the STB files indicated in
                % stbFilepathArray
                for i=m:-1:1
                    for j=n:-1:1
                        % Generate Particle, Track and TimeStep objects
                        % from the STB file
                        fprintf(['Importing STB data from file ',...
                            '[%d/%d]\n-> %s\n'],n-j+1+(n-j+1)*(m-i),m*n,...
                            stbFilepathArray{i,j})
                        [particleVector,trackVector,timeStepVector] =...
                            loadStbDataDavis10(stbFilepathArray{i,j});
                        % Assign generated objects to properties
                        obj(i,j).RootFile = stbFilepathArray{i,j};
                        obj(i,j).ParticleVector = particleVector;
                        obj(i,j).TimeStepVector = timeStepVector;
                        obj(i,j).TrackVector = trackVector;
                    end
                end
            end
        end
        
        %% Time average binning
        function [binningInputArray,binnedData] = timeAverageBinning(...
                obj,binningInputStruct)
            %timeAverageBinning Perform time averaged binning on the STB
            %dataset given by the current object or vector of objects
            %   [binningInputArray,binnedData,referenceCyclePeriod] =...
            %   bin(obj,binningInputStruct)
            % Assemble an array with all StbParticle objects contained in
            % the considered StbFlowRun objects
            particle4BinningArray = vertcat(obj.ParticleVector);
            % Perofrm binning of selected particle set
            [binnedData,binningInputArray] =...
                particle4BinningArray.bin(binningInputStruct);
        end
        
        %% Phase average binning
        function [binningInputArray,binnedData] = phaseAverageBinning(...
                obj,binningInputStruct,referencePeriod,timeOffset)
            %phaseAverageBinning Perform phase averaged binning on the STB
            %dataset given by the current object or vector of objects
            %   [binningInputArray,binnedData,referenceCyclePeriod] =...
            %   phaseAverageBinning(obj,binningInputStruct,...
            %   referencePeriod,timeOffset) where referencePeriod and
            %   timeOffset must have the same size of the StbRun object
            %   array considered. The first indicates the period of the
            %   average cycle of each StbRun object and is a mandatory
            %   input variable, the second indicates the temporal offset
            %   that has to be applied to the STB instants such that all
            %   runs start with the same phase.
            % Set phase average flag to true in binningInputStruct
            fprintf('Phase average binning: setting flag to true\n')
            binningInputStruct.phaseAverageFlag = true;
            % If time offset is not specified then assume zero time offset
            % for all STB runs
            if nargin<4
                timeOffset = zeros(size(obj));
            end
            % Find array with particles' positional and velocity info and
            % vector with particles' nondimensional time
            [binningInputArray,particleNondimensionalTimeVector] =...
                obj.getPhaseAverageBinningInput(referencePeriod,...
                timeOffset);
            % Save the nondimensional time vector into binningInputStruct
            binningInputStruct.particleNondimensionalTimeVector=...
                particleNondimensionalTimeVector;
            % Call binning function (all variables must be in SI units)
            binnedData = binning(binningInputArray,...
                binningInputStruct.binSize,binningInputStruct);
        end
        
        %% Get phase average binning input
        function [particleXyzUvwArray,...
                particleNondimensionalTimeVector] =...
                getPhaseAverageBinningInput(obj,referencePeriod,timeOffset)
            % Iterate through the StbRun objects to define the
            % nondimensional time instants of STB particles
            for i=length(obj):-1:1
                fprintf(['Setting reference period and time offset for',...
                    ' STB run [%d/%d]\n'],length(obj)-i+1,length(obj))
                % Assemble a vector with the measurement instants of
                % the particles of the current StbRun object,
                % zeroing such instants with the instant of the first
                % phase of the reference periodic signal
                particleTimeVector = vertcat(obj(i).ParticleVector.T)+...
                    timeOffset(i);
                % Nondimensionalize the time instants with the cycle
                % period and store result in a cell array where each
                % cell corresponds to a different StbRun object
                particleNondimensionalTimeArray{i,1} =...
                    particleTimeVector-floor(particleTimeVector/...
                    referencePeriod(i))*referencePeriod(i);
            end
            % Assemble a vector with the nondimensionalized time
            % instans from all the StbFlowRun considered
            particleNondimensionalTimeVector = cell2mat(...
                particleNondimensionalTimeArray);
            % Assemble an array with all StbParticle objects contained in
            % the considered StbFlowRun objects
            particle4BinningArray = vertcat(obj.ParticleVector);
            % Assemble the array with positional and velocity information
            particleXyzUvwArray =...
                particle4BinningArray.getBinningInputArray;
        end
        
        %% saveobj method
        function s = saveobj(obj)
            fprintf('Saving StbRun object\n')
            fprintf('-> path to original .dat file\n')
            s.rootFile = obj.RootFile;
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
            fprintf('-> array with particles'' time instant\n')
            s.timeVector = vertcat(obj.ParticleVector.T);
            fprintf('-> array with particles'' time step number\n')
            s.timeStepNoVector = vertcat(obj.ParticleVector.TimeStepNo);
            fprintf('-> array with particles'' track id\n')
            s.trackIdVector = vertcat(obj.ParticleVector.TrackId);
        end
    end
    
    %% loadobj method
    methods(Static)
        function obj = loadobj(s)
            newObj = StbRun;
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
            % Regenerate Track object
            fprintf('Regenerating array of StbTrack objects\n')
            trackVector = StbTrack(num2cell(uniqueTrackIdVector),...
                particlePerTrackArray);
            % Assign generated objects to properties
            newObj.RootFile = s.rootFile;
            newObj.ParticleVector = particleVector;
            newObj.TimeStepVector = timeStepObjVector;
            newObj.TrackVector = trackVector;
            obj = newObj;
        end
    end
end

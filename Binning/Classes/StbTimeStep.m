classdef StbTimeStep < matlab.mixin.Copyable
    %TimeStep Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        No
        T   % [s]
        ParticleVector
    end
    
    methods
        %% Constructor
        function obj = StbTimeStep(timeStepNoArray,instantArray,...
                particleArray)
            %StbTimeStep Construct an instance of this class
            %   Detailed explanation goes here
            % If number of input arguments is not zero then initialize the
            % object array with the size of the first input
            if nargin ~= 0
                [m,n] = size(timeStepNoArray);
                obj(m,n) = StbTimeStep;
                % Iterate through the elements of the input variables
                for i = m:-1:1
                    for j = n:-1:1
                        obj(i,j).No = timeStepNoArray(i,j);
                        obj(i,j).T = instantArray(i,j);
                        obj(i,j).ParticleVector = particleArray{i,j};
                    end
                end
            end
        end
        
        %% ParticleVector set method
        function set.ParticleVector(obj,particleVector)
            obj.ParticleVector = particleVector;
            % Assign current TimeStep object as parent of all StbParticle
            % objects in input
            for i = 1:length(obj.ParticleVector)
                obj.ParticleVector(i).TimeStep = obj;
            end
        end
        
        %% getTrackVector method
        function orderedTrackVector = getTrackVector(obj)
            % Returns the complete set of Track objects belonging to the
            % TimeStep object considered. The returned set is ordered
            % according to the track id
            particleVector = vertcat(obj.ParticleVector);
            trackVector = unique(vertcat(particleVector.Track));
            [~,I] = sort(vertcat(trackVector.TrackId));
            orderedTrackVector = trackVector(I);
        end
    end
end

classdef StbParticle < matlab.mixin.Copyable
    %StbParticle Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        X   % [mm]
        Y   % [mm]
        Z   % [mm]
        I
        U   % [m/s]
        V   % [m/s]
        W   % [m/s]
        VelocityMagnitude   % [m/s]
        Ax  % [m^2/s]
        Ay  % [m^2/s]
        Az  % [m^2/s]
        AccelerationMagnitude   % [m^2/s]
        Track
        TimeStep
    end
    properties (Dependent,SetAccess=private)
        TrackId
        T
        TimeStepNo
    end
    
    methods
        %% Constructor
        function obj = StbParticle(xArray,yArray,zArray,iArray,uArray,...
                vArray,wArray,velocityMagnitudeArray,axArray,ayArray,...
                azArray,accelerationMagnitudeArray)
            %StbParticle Construct an instance of this class
            %   Detailed explanation goes here
            % If number of input arguments is not zero then initialize the
            % object array with the size of the first input
            if nargin ~= 0
                [m,n] = size(xArray);
                obj(m,n) = StbParticle;
                % Iterate through the elements of the input variables
                for i = m:-1:1
                    for j = n:-1:1
                        obj(i,j).X = xArray(i,j);
                        obj(i,j).Y = yArray(i,j);
                        obj(i,j).Z = zArray(i,j);
                        obj(i,j).I = iArray(i,j);
                        obj(i,j).U = uArray(i,j);
                        obj(i,j).V = vArray(i,j);
                        obj(i,j).W = wArray(i,j);
                        obj(i,j).VelocityMagnitude =...
                            velocityMagnitudeArray(i,j);
                        obj(i,j).Ax = axArray(i,j);
                        obj(i,j).Ay = ayArray(i,j);
                        obj(i,j).Az = azArray(i,j);
                        obj(i,j).AccelerationMagnitude =...
                            accelerationMagnitudeArray(i,j);
                    end
                end
            end
        end
                
        %% TrackId, T and TimeStepNo get methods
        function trackId = get.TrackId(obj)
            trackId = obj.Track.TrackId;
        end
        function t = get.T(obj)
            t = obj.TimeStep.T;
        end
        function timeStepNo = get.TimeStepNo(obj)
            timeStepNo = obj.TimeStep.No;
        end
        
        %% plotParticle method
        function h = plot(obj,varargin)
            % Create an InputParser object
            p = inputParser;
            % Add inputs to the parsing scheme
            defaultColor = [0,.75,.75];
            addRequired(p,'obj',@(obj)isa(obj,'StbParticle'));
            addParameter(p,'color',defaultColor,@(x)...
                isnumeric(x)||ischar(x))
            addParameter(p,'targetAxes',gca)
            % Set properties to adjust parsing
            p.KeepUnmatched = true;
            % Parse the inputs
            parse(p,obj,varargin{:})
            % Generate scatter plot
            h = scatter3(p.Results.targetAxes,[obj.X],[obj.Y],[obj.Z],...
                'MarkerEdgeColor','k','MarkerFaceColor',p.Results.color);
            % Adjust plot
            makePlotNicer(struct('txtXlabel','$x$ [m]',...
                'txtYlabel','$y$ [m]','txtZlabel','$z$ [m]'))
            axis image
        end
        
        %% getBinningInputArray method
        function particleXyzUvwArray = getBinningInputArray(obj)
            % Assemble array with positional and velocity information of
            % all particles
            fprintf(['Generating array with positional and velocity ',...
                'information of particles\n'])
            particleXyzUvwArray = [vertcat(obj.X),vertcat(obj.Y),...
                vertcat(obj.Z),vertcat(obj.U),vertcat(obj.V),...
                vertcat(obj.W)];
        end
        
        %% bin method
        function [binnedData,particleXyzUvwArray] = bin(obj,...
                binningInputStruct)
            %bin Perform binning of the considered particle set.
            %   [binnedData,particleXyzUvwArray] = bin(obj,
            %   binningInputStruct) carries out an ensemble average of the 
            %   Lagrangian particle tracks given by the selected objects.
            %   binningInputStruct is a structure containing all kinds of
            %   input parameters for the ensemble average process (such as
            %   bin size, the overlap factor, the minimum number of
            %   particles, the averaging method, options for phase average
            %   and pressure calculation, etc.).
            % Assemble array with positional and velocity information of
            % all particles
            particleXyzUvwArray = obj.getBinningInputArray;
            % Call binning function  (all variables must be in SI units)
            binnedData = binning(particleXyzUvwArray,...
                binningInputStruct.binSize,binningInputStruct);
        end
    end
end

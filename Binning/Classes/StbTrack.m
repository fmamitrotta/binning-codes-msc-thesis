classdef StbTrack < matlab.mixin.Copyable
    %Track Summary of this class goes here
    %   Detailed explanation goes here
    
    %% Properties
    properties
        TrackId
    end
    properties (Dependent)
        T
        TimeStepNo
        X
        Y
        Z
        U
        V
        W
        ParticleVector
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
        ParticleVectorProxy
    end
    
    methods
        %% Constructor
        function obj = StbTrack(trackIdArray,particleArray)
            %StbTrack Construct an instance of this class
            %   obj = StbTrack(trackIdArray,particleArray) generates an
            %   array of StbTrack objects of the same size of trackIdArray.
            %   Both input are expected as cell arrays and they must have
            %   same size. trackIdArray contains the Id of the tracks and
            %   particleArray contains the StbParticle array corresponding
            %   to each track.
            % If number of input arguments is not zero then initialize the
            % object array with the size of the first input
            if nargin ~= 0
                [m,n] = size(trackIdArray);
                obj(m,n) = StbTrack;
                % Assign properties
                [obj.TrackId] = trackIdArray{:};
                [obj.ParticleVector] = particleArray{:};   
            end
        end
        
        %% ParticleVector set method
        function set.ParticleVector(obj,particleVector)
            % Assign ParticleVectory property
            obj.ParticleVectorProxy = particleVector;
            % Preallocate space for Tproxy and TimeStepNoProxy properties
            for i=length(particleVector):-1:1
                % Assign current Track object as parent of all StbParticle
                % objects in input
                particleVector(i).Track = obj;
                % Assign proxy properties
%                 obj.Tproxy(i) = particleVector(i).TimeStep.T;
%                 obj.TimeStepNoProxy(i) = particleVector(i).TimeStep.No;
%                 obj.Xproxy(i) = particleVector(i).X;
%                 obj.Yproxy(i) = particleVector(i).Y;
%                 obj.Zproxy(i) = particleVector(i).Z;
%                 obj.Uproxy(i) = particleVector(i).U;
%                 obj.Vproxy(i) = particleVector(i).V;
%                 obj.Wproxy(i) = particleVector(i).W;
            end
            timeStepVector = vertcat(particleVector.TimeStep);
            obj.Tproxy = vertcat(timeStepVector.T);
            obj.TimeStepNoProxy = vertcat(timeStepVector.No);
                obj.Xproxy = vertcat(particleVector.X);
                obj.Yproxy = vertcat(particleVector.Y);
                obj.Zproxy = vertcat(particleVector.Z);
                obj.Uproxy = vertcat(particleVector.U);
                obj.Vproxy = vertcat(particleVector.V);
                obj.Wproxy = vertcat(particleVector.W);
        end
        
        %% ParticleVector get method
        function particleVector = get.ParticleVector(obj)
            particleVector = obj.ParticleVectorProxy;
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
        
        %% plotTrack method
        function h = plot(obj,varargin)
            % Generate vector of StbParticle objects belonging to StbTrack
            % objects
            completeParticleVector = vertcat(obj.ParticleVector);
            % Call StbParticle plot method
            h = completeParticleVector.plot(varargin{:});
        end
    end
end


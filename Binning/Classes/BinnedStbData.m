classdef BinnedStbData < matlab.mixin.Copyable
    %BinnedStbData Eulerian description of a flow field obtained from the
    %binning of STB data.
    %   BinnedStbData contains the results of a binning operation carried
    %   out on a set of Lagrangian particle tracks, which form an Eulerian
    %   description of a flow field. A single object represents a
    %   time-averaged flow, while a vector of objects represents a
    %   phase-averaged flow, with each element corresponding to a different
    %   phase.
    
    properties
        NondimensionalTime  % double indicating the nondimensional time of the binned dataset, empty for time averaged results
        Xarray  % 3D double array with x coordinate of bins' centroid
        Yarray  % 3D double array with y coordinate of bins' centroid
        Zarray  % 3D double array with z coordinate of bins' centroid
        UmeanArray  % 3D double array with mean x velocity component
        VmeanArray  % 3D double array with mean y velocity component
        WmeanArray  % 3D double array with mean z velocity component
        UstdArray   % 3D double array with std of x velocity component
        VstdArray   % 3D double array with std of y velocity component
        WstdArray   % 3D double array with std of z velocity component
        NoParticlesArray    % 3D double array with number of particles in each bin
        % 3D double arrays containing Reynolds stresses in each bin
        UprimeUprimeArray
        VprimeVprimeArray
        WprimeWprimeArray
        UprimeVprimeArray
        UprimeWprimeArray
        VprimeWprimeArray
    end
    
    methods
        function obj = BinnedStbData(binCentroidXVector,...
                binCentroidYVector,binCentroidZVector,...
                velocityMeanArray,velocityStdArray,...
                noFoundParticlesArray,uprimeUprimeArray,...
                vprimeVprimeArray,wprimeWprimeArray,...
                uprimeVprimeArray,uprimeWprimeArray,...
                vprimeWprimeArray,binNondimensionalTimeVector)
            %BinnedStbData Construct an instance of this class
            %   obj = BinnedStbData(binCentroidXVector,binCentroidYVector,
            %   binCentroidZVector,velocityMeanArray,velocityStdArray,
            %   noFoundParticlesArray,uprimeUprimeArray,vprimeVprimeArray,
            %   wprimeWprimeArray,uprimeVprimeArray,uprimeWprimeArray,
            %   vprimeWprimeArray,binNondimensionalTimeVector) returns a
            %   BinnedStbData object or a vector of objects. If
            %   binNondimensionalTimeVector is given as input, the
            %   constructor expects all the other input variables to be
            %   cell array, with each cell corresponding to a different
            %   phase of the average cycle. Otherwise, all the input
            %   variables are expeceted as double arrays.
            % If number of input variables is different from zero proceed
            % with property assignment, otherwise generate empty object
            if nargin ~= 0
                % Generate the 3D arrays of the centroids' coordinates
                [xArray,yArray,zArray] = meshgrid(binCentroidXVector,...
                    binCentroidYVector,binCentroidZVector);
                % Distinguish between time averaged and phase averaged
                % results
                if nargin > 12
                    % If binNondimensionalTimeVector is given as input
                    % variable, then consider phase averaged results
                    % Generate an array of empty BinnedStbData objects of
                    % the same size of binNondimensionalTimeVector
                    [m,n] = size(binNondimensionalTimeVector);
                    obj(m,n) = BinnedStbData;
                    % Iterate through the elements of
                    % binNondimensionalTimeVector (the phases)
                    for i=1:numel(binNondimensionalTimeVector)
                        % Assign property to each element of the object
                        % array
                        obj(i).NondimensionalTime =...
                            binNondimensionalTimeVector(i);
                        obj(i).Xarray = xArray;
                        obj(i).Yarray = yArray;
                        obj(i).Zarray = zArray;
                        obj(i).UmeanArray = velocityMeanArray{i}(:,:,:,1);
                        obj(i).VmeanArray = velocityMeanArray{i}(:,:,:,2);
                        obj(i).WmeanArray = velocityMeanArray{i}(:,:,:,3);
                        obj(i).UstdArray = velocityStdArray{i}(:,:,:,1);
                        obj(i).VstdArray = velocityStdArray{i}(:,:,:,2);
                        obj(i).WstdArray = velocityStdArray{i}(:,:,:,3);
                        obj(i).NoParticlesArray = noFoundParticlesArray{i};
                        obj(i).UprimeUprimeArray = uprimeUprimeArray{i};
                        obj(i).VprimeVprimeArray = vprimeVprimeArray{i};
                        obj(i).WprimeWprimeArray = wprimeWprimeArray{i};
                        obj(i).UprimeVprimeArray = uprimeVprimeArray{i};
                        obj(i).UprimeWprimeArray = uprimeWprimeArray{i};
                        obj(i).VprimeWprimeArray = vprimeWprimeArray{i};
                    end
                else
                    % If binNondimensionalTimeVector is not given as input
                    % variable, then consider time averaged results (only
                    % one object is generated)
                    % Assign properties to object
                    obj.Xarray = xArray;
                    obj.Yarray = yArray;
                    obj.Zarray = zArray;
                    obj.UmeanArray = velocityMeanArray(:,:,:,1);
                    obj.VmeanArray = velocityMeanArray(:,:,:,2);
                    obj.WmeanArray = velocityMeanArray(:,:,:,3);
                    obj.UstdArray = velocityStdArray(:,:,:,1);
                    obj.VstdArray = velocityStdArray(:,:,:,2);
                    obj.WstdArray = velocityStdArray(:,:,:,3);
                    obj.NoParticlesArray = noFoundParticlesArray;
                    obj.UprimeUprimeArray = uprimeUprimeArray;
                    obj.VprimeVprimeArray = vprimeVprimeArray;
                    obj.WprimeWprimeArray = wprimeWprimeArray;
                    obj.UprimeVprimeArray = uprimeVprimeArray;
                    obj.UprimeWprimeArray = uprimeWprimeArray;
                    obj.VprimeWprimeArray = vprimeWprimeArray;
                end
            end
        end
        
        %% Write binned data into a binary file for tecplot
        function writeTecplotBinary(obj,filename,folderPath)
            % Set default file name if not given as input
            if nargin==1
                filename = 'BinnedTracks.plt';
            elseif nargin>=2
                filename = [matlab.lang.makeValidName(filename),'.plt'];
                if nargin>=3
                    filename = fullfile(folderPath,filename);
                end
            end
            % Check input object and define a phase vector if lentgh>1
            if length(obj)>1
                nondimensionalTimeVector = [obj.NondimensionalTime];
            end
            % Initialize input structure for mat2tecplot
            % Set number of variables
            tdata.Nvar=17;
            % Set double precision for each variable
            tdata.vformat = 2*ones(1,tdata.Nvar);
            % Set variable names
            tdata.varnames={'X','Y','Z',...                             % 3
                'u [m/s]','v  [m/s]','w  [m/s]','|u|  [m/s]',...        % 4 (1-4)
                'u_std [m/s]','v_std  [m/s]','w_std  [m/s]',...         % 3 (5-7)
                'Om_x [1/s]','Om_y [1/s]','Om_z [1/s]','|Om| [1/s]',... % 4 (8-11)
                'Q','np','isvalid'};                                    % 3 (12-14)
            % Iterate through the ojbect vector
            for i=length(obj):-1:1
                % Add time related parameters if object represent a phase
                % averaged flow
                if length(obj)>1
                    % Set zone name with the nondimensional time instant
                    tdata.cubes(i).zonename = sprintf('t/Tref=%.2f',...
                        nondimensionalTimeVector(i));
                    % Set solution time of the zone
                    tdata.cubes(i).solutiontime =...
                        nondimensionalTimeVector(i);
                    % Set a unique Strand ID to group together zones
                    % belonging to different time steps
                    tdata.cubes(i).strandID = uint8(1);
                end
                % Extract vector containing all velocities in the grid
                uVector = obj(i).UmeanArray(:);
                vVector = obj(i).VmeanArray(:);
                wVector = obj(i).WmeanArray(:);
                % Generate 3D array with velocity magnitude
                velocityMagnitudeArray = (obj(i).UmeanArray.^2+...
                    obj(i).VmeanArray.^2+obj(i).WmeanArray.^2).^(0.5);
                % Determine the vector spacing in all directions
                xVectorSpacing = obj(i).Xarray(1,2,1)-obj(i).Xarray(1,1,1);
                yVectorSpacing = obj(i).Yarray(2,1,1)-obj(i).Yarray(1,1,1);
                zVectorSpacing = obj(i).Zarray(1,1,2)-obj(i).Zarray(1,1,1);
                % Assemble a 4D array containing all velocities components
                % in the grid
                overallVelocityArray(:,:,:,3) = obj(i).WmeanArray;
                overallVelocityArray(:,:,:,2) = obj(i).VmeanArray;
                overallVelocityArray(:,:,:,1) = obj(i).UmeanArray;
                % Calculate vorticity field and Q criterion
                omegaArray = xicalc3d(xVectorSpacing,yVectorSpacing,...
                    zVectorSpacing,overallVelocityArray,200);
                qCriterionArray  = qcritcalc(xVectorSpacing,...
                    yVectorSpacing,zVectorSpacing,overallVelocityArray);
                % Clear overallVelocityArray for next iteration
                clear overallVelocityArray
                % Retrieve the components of the vorticity field
                omegaXarray = omegaArray(:,:,:,1);
                omegaYarray = omegaArray(:,:,:,2);
                omegaZarray = omegaArray(:,:,:,3);
                % Calculate the magnitude of the vorticity field
                omegaMagnitudeArray = sqrt(sum(omegaArray.^2,4));
                % Find all NaN velocities
                nanmask = isnan(uVector) | isnan(vVector) | isnan(wVector);
                % Generate a 3D array to flag invalid cells in the grid
                isValid = reshape(~nanmask,size(obj(i).Xarray));
                % Grid coordinates
                tdata.cubes(i).x = obj(i).Xarray; % x
                tdata.cubes(i).y = obj(i).Yarray; % y
                tdata.cubes(i).z = obj(i).Zarray; % z
                % Mean velocity
                tdata.cubes(i).v(1,:,:,: ) = obj(i).UmeanArray; % u
                tdata.cubes(i).v(2,:,:,: ) = obj(i).VmeanArray; % v
                tdata.cubes(i).v(3,:,:,: ) = obj(i).WmeanArray; % w
                tdata.cubes(i).v(4,:,:,: ) = velocityMagnitudeArray;   %|u|
                % Velocity standard deviation
                tdata.cubes(i).v(5,:,:,: ) = obj(i).UstdArray; % u_std
                tdata.cubes(i).v(6,:,:,: ) = obj(i).VstdArray; % v_std
                tdata.cubes(i).v(7,:,:,: ) = obj(i).WstdArray; % w_std
                % Vorticity
                tdata.cubes(i).v(8,:,:,: ) = omegaXarray;
                tdata.cubes(i).v(9,:,:,: ) = omegaYarray;
                tdata.cubes(i).v(10,:,:,: ) = omegaZarray;
                tdata.cubes(i).v(11,:,:,: ) = omegaMagnitudeArray;
                % Q-criterion
                tdata.cubes(i).v(12,:,:,: ) = qCriterionArray;
                % Number of particles
                tdata.cubes(i).v(13,:,:,: ) = obj(i).NoParticlesArray;
                % NaN flag
                tdata.cubes(i).v(14,:,:,: ) = isValid;
            end
            % Write tecplot file
            mat2tecplot(tdata,filename);
        end
    end
end

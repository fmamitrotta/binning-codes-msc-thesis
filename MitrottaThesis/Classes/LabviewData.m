classdef LabviewData < matlab.mixin.Copyable
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        RootFile = '';
        BatchNo
        TimeVector
        PivTriggerVector
        GustVane1Vector
        GustVane2Vector
        ActuationVector
        VibrometerVector
        GustTriggerVector
    end
    
    properties (Constant)
        PivTriggerThreshold = 4;  % [V]
        GustTriggerThreshold = 2;   % [V]
    end
    
    methods
        function obj = LabviewData(labviewFilepathArray)
            %UNTITLED Construct an instance of this class
            %   Detailed explanation goes here
            if nargin~=0
                switch class(labviewFilepathArray)
                    case 'cell'
                        % If stbFilepathArray is a cell array initialize
                        % object with same dimension as stbFilepathArray
                        [m,n] = size(labviewFilepathArray);
                    case 'char'
                        % If stbFilepathArray is a char initialize object
                        % as 1x1 array
                        m = 1;
                        n = 1;
                        labviewFilepathArray =...
                            textscan(labviewFilepathArray,'%s',...
                            'Delimiter','\n');
                        labviewFilepathArray = labviewFilepathArray{:};
                    otherwise
                        % Throw error if stbFilepathArray is neither a cell
                        % nor a char
                        error(['Error. \nInput must be a cell or a ',...
                            'char, not a %s.'],class(labviewFilepathArray))
                end
                % Initialize a temporary cell array
                tempCellArray = cell(m,n);
                % Iterate through the STB files indicated in
                % stbFilepathArray
                for i=m:-1:1
                    for j=n:-1:1
                        % Load data from file
                        fprintf(['Importing LabVIEW data from file ',...
                            '[%d/%d]\n-> %s\n'],n-j+1+(n-j+1)*(m-i),m*n,...
                            labviewFilepathArray{i,j})
                        A = importdata(labviewFilepathArray{i,j});
                        % Columns in A represent the following
                        % 1 Time
                        % 2 Camera_trig
                        % 3 Laser1
                        % 4 Laser2
                        % 5 Actuation
                        % 6 Vibrometer_trig
                        % 7 Gust_trig
                        % Find number of data batches
                        timeDiffVector = diff(A(:,1));
                        timeJumpIndexVector =...
                            find(timeDiffVector>2*mean(timeDiffVector));
                        if ~isempty(timeJumpIndexVector)
                            noDataBatches = length(timeJumpIndexVector)+1;
                        else
                            noDataBatches = 1;
                        end
                        fprintf('-> %d data batches found\n',noDataBatches)
                        % Find indices of batch start and end
                        batchStartIndexVector = [1;timeJumpIndexVector+1];
                        batchEndIndexVector =...
                            [timeJumpIndexVector;size(A,1)];
                        % Generate an object vector with noDataBatches
                        % elements
                        tempCellArray{i,j}(noDataBatches,1) = LabviewData;
                        % Assign root filepath
                        C = cell(noDataBatches,1);
                        C(:) = labviewFilepathArray(i,j);
                        [tempCellArray{i,j}.RootFile] = C{:};
                        % Assign batch number
                        C = num2cell(1:noDataBatches);
                        [tempCellArray{i,j}.BatchNo] = C{:};
                        % Assign time vector
                        timeArray = arrayfun(@(x) A(...
                            batchStartIndexVector(x):...
                            batchEndIndexVector(x),1),...
                            1:noDataBatches,'UniformOutput',false);
                        [tempCellArray{i,j}.TimeVector] = timeArray{:};
                        % Assign piv trigger vector
                        pivTriggerArray = arrayfun(@(x) A(...
                            batchStartIndexVector(x):...
                            batchEndIndexVector(x),2),...
                            1:noDataBatches,'UniformOutput',false);
                        [tempCellArray{i,j}.PivTriggerVector] =...
                            pivTriggerArray{:};
                        % Assign gust vane 1 vector
                        gustVane1Array = arrayfun(@(x)...
                            calibrationGustVane1(A(...
                            batchStartIndexVector(x):...
                            batchEndIndexVector(x),3)),...
                            1:noDataBatches,'UniformOutput',false);
                        [tempCellArray{i,j}.GustVane1Vector] =...
                            gustVane1Array{:};
                        % Assign gust vane 2 vector
                        gustVane2Array = arrayfun(@(x)...
                            calibrationGustVane2(A(...
                            batchStartIndexVector(x):...
                            batchEndIndexVector(x),4)),...
                            1:noDataBatches,'UniformOutput',false);
                        [tempCellArray{i,j}.GustVane2Vector] =...
                            gustVane2Array{:};
                        % Assign actuation vector
                        actuationArray = arrayfun(@(x) A(...
                            batchStartIndexVector(x):...
                            batchEndIndexVector(x),5),...
                            1:noDataBatches,'UniformOutput',false);
                        [tempCellArray{i,j}.ActuationVector] =...
                            actuationArray{:};
                        % Assign vibrometer vector
                        vibrometerArray = arrayfun(@(x) A(...
                            batchStartIndexVector(x):...
                            batchEndIndexVector(x),6),...
                            1:noDataBatches,'UniformOutput',false);
                        [tempCellArray{i,j}.VibrometerVector] =...
                            vibrometerArray{:};
                        % Assign gust trigger vector
                        gustTriggerArray = arrayfun(@(x) A(...
                            batchStartIndexVector(x):...
                            batchEndIndexVector(x),7),...
                            1:noDataBatches,'UniformOutput',false);
                        [tempCellArray{i,j}.GustTriggerVector] =...
                            gustTriggerArray{:};
                    end
                end
                % Assemble output object
                obj = vertcat(tempCellArray{:});
            end
        end
        
        %% Find PIV instants
        function pivTimeVector = findPivAcquisitionInstants(obj)
            %findPivAcquisitionInstants Returns the vector of PIV image
            %acquisition instants.
            %   Detailed explanation goes here
            % Find all instants in the LabviewData object such that the
            % PivTriggerVector property is larger than the
            % PivTriggerThreshold property
            pivAboveThresholdTimeVector = obj.TimeVector(...
                obj.PivTriggerVector>=obj.PivTriggerThreshold);
            % Take the difference between the time instants found
            pivAboveThresholdTimeDiffVector =...
                diff(pivAboveThresholdTimeVector);
            % Find the indices of the time instants vector that separate
            % PIV trigger pulses
            triggerPulseSeparationIndexVector = find(...
                pivAboveThresholdTimeDiffVector>...
                mean(pivAboveThresholdTimeDiffVector)+...
                std(pivAboveThresholdTimeDiffVector));
            % Complement index vector with start and end values
            triggerPulseSeparationIndexVector = [0;...
                triggerPulseSeparationIndexVector;...
                length(pivAboveThresholdTimeDiffVector)];
            % Find PIV instants as the average time between the two time
            % instants defining each PIV trigger pulse
            pivTimeVector = mean(cell2mat(arrayfun(@(x)...
                [pivAboveThresholdTimeVector(...
                triggerPulseSeparationIndexVector(x)+1),...
                pivAboveThresholdTimeVector(...
                triggerPulseSeparationIndexVector(x+1))],...
                (1:length(triggerPulseSeparationIndexVector)-1)',...
                'UniformOutput',false)),2);
            % Return only the PIV instants coming after the gust trigger
            pivTimeVector = pivTimeVector(pivTimeVector>=obj.TimeVector(...
                find(obj.GustTriggerVector>=obj.GustTriggerThreshold,1)));
            %
            % Plot PIV acquisition instants (debugging)
            figure
            hold on
            plot(obj.TimeVector,obj.PivTriggerVector)
            scatter(obj.TimeVector(find(obj.GustTriggerVector>=...
                obj.GustTriggerThreshold,1)),obj.GustTriggerThreshold)
            scatter(pivTimeVector,ones(size(pivTimeVector))*...
                obj.PivTriggerThreshold,'g')
            %
        end
        
        %% Plot gust vanes angle
        function plotGustVanesAngle(obj,gustGeneratorFilePath)
            % If path to gust generator file is present, then correct the
            % mean value of the LabVIEW gust generator signals
            if nargin>1 && ~isempty(gustGeneratorFilePath)
                % Load gust generator data
                A = importdata(gustGeneratorFilePath);
                % Change decimal separator from comma to dot
                A = strrep(A,',','.');
                % Obtain a numerical array with the gust vane signals from
                % the gust generator data (columns 5 and 6)
                dataArray = cell2mat(cellfun(@(x) sscanf(x,...
                    [repmat('%f',1,6),repmat('%*f',1,13)])',...
                    A(3:end),'UniformOutput',false));
                % Calculate the correction term as the the mean of the
                % signals coming from the gust generator data minus the
                % mean from the gust vane signals acquired in LabVIEW
                gustVane1Correction = mean(dataArray(:,5))-...
                    mean(obj.GustVane1Vector);
                gustVane2Correction = mean(dataArray(:,6))-...
                    mean(obj.GustVane2Vector);
                % Correct mean of gust vanes signal
                gustVane1PlotVector = obj.GustVane1Vector+...
                    gustVane1Correction;
                gustVane2PlotVector = obj.GustVane2Vector+...
                    gustVane2Correction;
            else
                % Set correction terms to 0
                gustVane1PlotVector = obj.GustVane1Vector;
                gustVane2PlotVector = obj.GustVane2Vector;
            end
            % Plot gust vanes angle
            figure
            plot(obj.TimeVector-obj.TimeVector(1),gustVane1PlotVector,...
                obj.TimeVector-obj.TimeVector(1),gustVane2PlotVector)
            % Find the acquisition instants of PIV measurement (last
            % instant used to set the x-axis limits)
            pivTimeVector = obj.findPivAcquisitionInstants;
            % Make plot nicer
            specificationStruct = struct('txtXlabel','Time [s]',...
                'txtYlabel','Gust vanes angle [deg]',...
                'xlim',[0,pivTimeVector(end)-obj.TimeVector(1)],...
                'legendArray',{{'Gust vane 1','Gust vane 2'}});
            makePlotNicer(specificationStruct)
        end
        
        %% Find phase average reference
        function [phase0ReferenceTime,referencePeriod,noCycles,...
                pivTimeVector] = calculatePhaseAverageReference(obj,...
                expectedFrequency,gustGeneratorFilePath,plotFlag)
            % Find the acquisition instants of PIV measurement
            pivTimeVector = obj.findPivAcquisitionInstants;
            % Determine the time vector to be used in the sine fitting
            timeVector4SineFit = obj.TimeVector(...
                obj.TimeVector>=pivTimeVector(1)&...
                obj.TimeVector<=pivTimeVector(end));
            % Determine the part of the signals from the gust vanes to
            % be used in the sine fitting
            gustVane1Vector4SineFit =...
                obj.GustVane1Vector(...
                obj.TimeVector>=pivTimeVector(1)&...
                obj.TimeVector<=pivTimeVector(end));
            gustVane2Vector4SineFit =...
                obj.GustVane2Vector(...
                obj.TimeVector>=pivTimeVector(1)&...
                obj.TimeVector<=pivTimeVector(end));
            % If path to gust generator file is present, then correct the
            % mean value of the LabVIEW gust generator signals
            if nargin>2 && ~isempty(gustGeneratorFilePath)
                % Load gust generator data
                A = importdata(gustGeneratorFilePath);
                % Change decimal separator from comma to dot
                A = strrep(A,',','.');
                % Obtain a numerical array with the gust vane signals from
                % the gust generator data (columns 5 and 6)
                dataArray = cell2mat(cellfun(@(x) sscanf(x,...
                    [repmat('%f',1,6),repmat('%*f',1,13)])',...
                    A(3:end),'UniformOutput',false));
                % Calculate the correction term as the the mean of the
                % signals coming from the gust generator data minus the
                % mean from the gust vane signals acquired in LabVIEW
                gustVane1Correction = mean(dataArray(:,5))-...
                    mean(gustVane1Vector4SineFit);
                gustVane2Correction = mean(dataArray(:,6))-...
                    mean(gustVane2Vector4SineFit);
                % Correct mean of gust vanes signal
                gustVane1Vector4SineFit = gustVane1Vector4SineFit+...
                    gustVane1Correction;
                gustVane2Vector4SineFit = gustVane2Vector4SineFit+...
                    gustVane2Correction;
            else
                % Set correction terms to 0
                gustVane1Correction = 0;
                gustVane2Correction = 0;
            end
            % Calculate the maximum expected amplitude as the maximum
            % difference among the recorded points of the gust vanes'
            % signal (the average between the two gust vanes is
            % considered)
            maxExpectedAmplitude = mean([max(...
                gustVane1Vector4SineFit)-...
                min(gustVane1Vector4SineFit),...
                max(gustVane2Vector4SineFit)-...
                min(gustVane2Vector4SineFit)])/2;
            startOffset = mean([gustVane1Vector4SineFit;...
                gustVane2Vector4SineFit]);
            % Set method, coefficient boundaries and start point for
            % the fit
            maxOmega = 2*pi/(pivTimeVector(2)-pivTimeVector(1))/2;
            fo = fitoptions('Display','iter',...
                'Method','NonlinearLeastSquares',...
                'Lower',[maxExpectedAmplitude*0.5,0,0,...
                -2*abs(startOffset)],...
                'Upper',[maxExpectedAmplitude,maxOmega,2*pi,...
                2*abs(startOffset)],...
                'StartPoint',[maxExpectedAmplitude,...
                2*pi*expectedFrequency,0,startOffset]);
            % Generate fittype object using a sine wave
            ft = fittype('a*sin(b*x+c)+d',...
                'dependent',{'y'},'independent',{'x'},...
                'coefficients',{'a','b','c','d'},'options',fo);
            % Perform fitting on the gust vanes' signal
            fprintf('Fitting sine wave to signal of gust vanes\n')
            gustVanesFit = fit([timeVector4SineFit;timeVector4SineFit],...
                [gustVane1Vector4SineFit;gustVane2Vector4SineFit],ft);
            %
            % Plot fit result (debugging)
%             figure
            hold on
            plot(gustVanesFit,timeVector4SineFit,gustVane1Vector4SineFit)
            plot(gustVanesFit,timeVector4SineFit,gustVane2Vector4SineFit)
            %
            % Determine coefficients of the fitting
            coeffvalsGustVanes = coeffvalues(gustVanesFit);
            % Determine reference angular frequency, period and phase shift
            referenceOmega = coeffvalsGustVanes(2);
            referencePeriod = 2*pi/coeffvalsGustVanes(2);
            referencePhaseShift = coeffvalsGustVanes(3);
            % Calculate the reference time for phase 0 taking the closest
            % multiple of 2*pi to the first piv recording instant
            phase0ReferenceTime = (floor((referenceOmega*...
                pivTimeVector(1)+referencePhaseShift)/(2*pi))*2*pi-...
                referencePhaseShift)/referenceOmega;
            %
            % Plot phase 0 reference time (debugging)
            plot(ones(1,2)*phase0ReferenceTime,...
                [-obj.PivTriggerThreshold,obj.PivTriggerThreshold])
            %
            % Calculate the number of cycles included in the current
            % StbRun object
            noCycles = ceil((pivTimeVector(end)-...
                phase0ReferenceTime)/referencePeriod);
            % If flag is set, plot fit with identified cycles
            if nargin>3 && plotFlag
                % Retrieve lines colormap array
                c = lines;
                % Generate subplots and axes for the two gust vanes
                figure
                [ha,~] = tight_subplot(2,1,.05,.15,.12);
                gustVane1Axes = ha(1);
                hold(gustVane1Axes,'on')
                gustVane2Axes = ha(2);
                hold(gustVane2Axes,'on')
                % Set current axes as the one of gust vane 1
                fig = gcf;
                fig.CurrentAxes = gustVane1Axes;
                % Plot the fit and the data
                hFit = plot(gustVanesFit,...
                    obj.TimeVector-phase0ReferenceTime,...
                    obj.GustVane1Vector+gustVane1Correction);
                % Retrieve y-axis limits
                yl1 = ylim(gustVane1Axes);
                % Plot first and last PIV instants
                hPivInstants = plot([ones(2,1)*pivTimeVector(1),...
                    ones(2,1)*pivTimeVector(end)]-phase0ReferenceTime,...
                    yl1'*ones(1,2),'--','color',c(5,:));
                % Plot cycles delimitation
                hCycles = plot(cell2mat(arrayfun(@(x)...
                    phase0ReferenceTime+referencePeriod*x*ones(2,1),...
                    0:noCycles,'UniformOutput',false))-...
                    phase0ReferenceTime,...
                    repmat([yl1(1);yl1(2)],1,noCycles+1),'--',...
                    'color',c(2,:));
                % Set current axes as the one of gust vane 2
                fig.CurrentAxes = gustVane2Axes;
                % Plot the fit and the data
                plot(gustVanesFit,obj.TimeVector-phase0ReferenceTime,...
                    obj.GustVane2Vector+gustVane2Correction)
                % Retrieve y-axis limits
                yl2 = ylim(gustVane2Axes);
                % Plot first and last PIV instants
                plot([ones(2,1)*pivTimeVector(1),...
                    ones(2,1)*pivTimeVector(end)]-phase0ReferenceTime,...
                    yl2'*ones(1,2),'--','color',c(5,:));
                % Plot cycles delimitation
                plot(cell2mat(arrayfun(@(x)...
                    phase0ReferenceTime+referencePeriod*x*ones(2,1),...
                    0:noCycles,'UniformOutput',false))-...
                    phase0ReferenceTime,...
                    repmat([yl2(1);yl2(2)],1,noCycles+1),'--',...
                    'color',c(2,:))
                % Set axes of subplots
                set(gustVane1Axes,'XLimSpec','tight')
                set(gustVane2Axes,'XLimSpec','tight')
                linkaxes([gustVane1Axes,gustVane2Axes],'x')
                xticklabels(gustVane2Axes,'auto')
                xticklabels(gustVane1Axes,{})
                yticklabels(gustVane2Axes,'auto')
                yticklabels(gustVane1Axes,'auto')
                set(gustVane1Axes,'xlabel',[])
                % Make gust vane 1 plot nicer
                plotSpecificationStruct = struct(...
                    'targetAxes',gustVane1Axes,...
                    'txtYlabel',{{'Gust vane 1';'angle [deg]'}},...
                    'ylim',yl1,...
                    'lineHandleVector',[hFit;hPivInstants(1);hCycles(1)],...
                    'legendArray',{{'Measured data','Fitted sine wave',...
                    'First and last PIV images','Cycles'' bounds'}});
                makePlotNicer(plotSpecificationStruct)
                % Make gust vane 2 plot nicer
                plotSpecificationStruct = struct(...
                    'targetAxes',gustVane2Axes,...
                    'txtXlabel','Time [s]',...
                    'txtYlabel',{{'Gust vane 2';'angle [deg]'}},...
                    'ylim',yl2);
                makePlotNicer(plotSpecificationStruct)
                legend(gustVane2Axes,'off')
            end
        end
    end
end

function angleVector = calibrationGustVane1(voltageVector)
voltage2FitVector = [3.4635
    3.4520
    3.4489
    3.4410
    3.4398
    3.4193
    3.4205
    3.4010
    3.4042
    3.3871
    3.3862
    3.3751
    3.3754
    3.3589
    3.3493
    3.3480
    3.3355
    3.3212
    3.3197
    3.3153
    3.3048
    3.3015
    3.2937
    3.2848
    3.2763
    3.2617
    3.2561
    3.2542
    3.2459
    3.2353
    3.2281
    3.2268
    3.2212
    3.2091
    3.2039
    3.1935
    3.1868
    3.1794
    3.1753
    3.1655
    3.1549
    3.1539
    3.1312
    3.1254
    3.1176
    3.1135
    3.1041
    3.1020
    3.0844
    3.0798
    3.0664
    3.0580
    3.0538
    3.0472
    3.0364
    3.0318
    3.0165
    3.0135
    3.0043
    2.9984
    2.9967
    2.9898
    2.9745
    2.9614
    2.9550
    2.9445
    2.9409
    2.9364
    2.9243
    2.9164
    2.9018
    2.9003
    2.8979
    2.8886
    2.8815
    2.8602
    2.8618
    2.8489
    2.8430
    2.8335
    2.8230
    2.8190
    2.8103
    2.7938
    2.7856
    2.7802
    2.7732
    2.7691
    2.7643
    2.7529
    2.7424
    2.7306
    2.7239
    2.7168
    2.7091
    2.6991
    2.6903
    2.6919
    2.6985
    2.7097
    2.7170
    2.7237
    2.7329
    2.7421
    2.7541
    2.7648
    2.7687
    2.7727
    2.7833
    2.7947
    2.8033
    2.8102
    2.8187
    2.8282
    2.8323
    2.8445
    2.8475
    2.8616
    2.8690
    2.8811
    2.8865
    2.8968
    2.9000
    2.9169
    2.9204
    2.9262
    2.9376
    2.9414
    2.9483
    2.9533
    2.9608
    2.9630
    2.9794
    2.9889
    2.9960
    3.0014
    3.0068
    3.0137
    3.0249
    3.0338
    3.0422
    3.0544
    3.0584
    3.0651
    3.0745
    3.0849
    3.0966
    3.0985
    3.1125
    3.1135
    3.1283
    3.1329
    3.1357
    3.1471
    3.1483
    3.1654
    3.1660
    3.1845
    3.1861
    3.1907
    3.2035
    3.2072
    3.2137
    3.2265
    3.2310
    3.2459
    3.2475
    3.2593
    3.2651
    3.2766
    3.2850
    3.2938
    3.2939
    3.3051
    3.3132
    3.3233
    3.3263
    3.3351
    3.3454
    3.3541
    3.3599
    3.3688
    3.3751
    3.3848
    3.3215
    3.2349
    3.4054
    3.4080
    3.4259
    3.4355
    3.4355
    3.4486
    3.4556
    3.4661];
angle2FitVector = [-12.0006
    -11.7492
    -11.4996
    -11.2500
    -11.0004
    -10.7508
    -10.4994
    -10.2498
    -10.0002
    -9.7506
    -9.4992
    -9.2496
    -9.0000
    -8.7504
    -8.5008
    -8.2494
    -7.9998
    -7.7502
    -7.5006
    -7.2492
    -6.9996
    -6.7500
    -6.5004
    -6.2508
    -5.9994
    -5.7498
    -5.5002
    -5.2506
    -4.9992
    -4.7496
    -4.5000
    -4.2504
    -4.0008
    -3.7494
    -3.4998
    -3.2502
    -3.0006
    -2.7492
    -2.4996
    -2.2500
    -2.0004
    0.0000
    -2.0004
    -1.4994
    -1.2498
    -1.0002
    -0.7506
    -0.4992
    -0.2496
    0.0000
    0.2496
    0.4992
    0.7506
    1.0002
    1.2498
    1.4994
    1.7508
    2.0004
    2.2500
    2.4996
    2.7492
    3.0006
    3.2502
    3.4998
    3.7494
    4.0008
    4.2504
    4.5000
    4.7496
    4.9992
    5.2506
    5.5002
    5.7498
    5.9994
    6.2508
    6.5004
    6.7500
    6.9996
    7.2492
    7.5006
    7.7502
    7.9998
    8.2494
    8.5008
    8.7504
    9.0000
    9.2496
    9.4992
    9.7506
    10.0002
    10.2498
    10.4994
    10.7508
    11.0004
    11.2500
    11.4996
    11.7492
    12.0006
    12.0006
    11.7492
    11.4996
    11.2500
    11.0004
    10.7508
    10.4994
    10.2498
    10.0002
    9.7506
    9.4992
    9.2496
    9.0000
    8.7504
    8.5008
    8.2494
    7.9998
    7.7502
    7.5006
    7.2492
    6.9996
    6.7500
    6.5004
    6.2508
    5.9994
    5.7498
    5.5002
    5.2506
    4.9992
    4.7496
    4.5000
    4.2504
    4.0008
    3.7494
    3.4998
    3.2502
    3.0006
    2.7492
    2.4996
    2.2500
    2.0004
    1.7508
    1.4994
    1.2498
    1.0002
    0.7506
    0.4992
    0.2496
    0.0000
    -0.2496
    -0.4992
    -0.7506
    -1.0002
    -1.2498
    -1.4994
    -1.7508
    -2.0004
    -2.2500
    -2.4996
    -2.7492
    -3.0006
    -3.2502
    -3.4998
    -3.7494
    -4.0008
    -4.2504
    -4.5000
    -4.7496
    -4.9992
    -5.2506
    -5.5002
    -5.7498
    -5.9994
    -6.2508
    -6.5004
    -6.7500
    -6.9996
    -7.2492
    -7.5006
    -7.7502
    -7.9998
    -8.2494
    -8.5008
    -8.7504
    -9.0000
    -9.2496
    -9.4992
    -0.9615
    -10.0002
    -10.2498
    -10.4994
    -10.7508
    -11.0004
    -11.2500
    -11.4996
    -11.7492];
% Perform linear fitting
p = polyfit(voltage2FitVector,angle2FitVector,1);
% Evaluate fitting with input voltage vector
angleVector = polyval(p,voltageVector);
end

function angleVector = calibrationGustVane2(voltageVector)
voltage2FitVector = [2.7909
    2.7778
    2.7964
    2.7878
    2.8030
    2.8333
    2.8420
    2.8415
    2.8558
    2.8436
    2.8754
    2.8760
    2.8751
    2.8771
    2.8919
    2.8913
    2.9094
    2.9109
    2.9458
    2.9454
    2.9537
    2.9481
    2.9629
    2.9622
    2.9924
    2.9911
    3.0017
    3.0060
    3.0213
    3.0174
    3.0475
    3.0482
    3.0613
    3.0651
    3.0664
    3.0625
    3.0867
    3.0940
    3.1129
    3.1244
    3.1292
    3.1304
    3.1382
    3.1425
    3.1505
    3.1671
    3.1755
    3.1816
    3.1864
    3.1929
    3.2052
    3.2126
    3.2204
    3.2343
    3.2430
    3.2505
    3.2593
    3.2616
    3.2679
    3.2765
    3.2855
    3.2969
    3.3102
    3.3240
    3.3295
    3.3303
    3.3387
    3.3443
    3.3518
    3.3683
    3.3770
    3.3904
    3.4022
    3.4007
    3.4084
    3.4147
    3.4257
    3.4328
    3.4416
    3.4549
    3.4639
    3.4711
    3.4725
    3.4824
    3.4908
    3.4991
    3.5094
    3.5194
    3.5285
    3.5358
    3.5417
    3.5502
    3.5573
    3.5666
    3.5704
    3.5780
    3.5857
    3.5904
    3.5796
    3.5720
    3.5703
    3.5682
    3.5556
    3.5553
    3.5412
    3.5395
    3.5335
    3.5258
    3.5135
    3.5068
    3.4918
    3.4923
    3.4761
    3.4721
    3.4642
    3.4646
    3.4515
    3.4408
    3.4345
    3.4163
    3.4136
    3.4061
    3.3985
    3.3989
    3.3895
    3.3732
    3.3684
    3.3513
    3.3441
    3.3394
    3.3300
    3.3314
    3.3245
    3.3100
    3.2962
    3.2864
    3.2755
    3.2676
    3.2635
    3.2605
    3.2511
    3.2462
    3.2348
    3.2188
    3.2051
    3.2028
    3.1968
    3.1888
    3.1904
    3.1745
    3.1695
    3.1504
    3.1434
    3.1329
    3.1305
    3.1315
    3.1244
    3.1111
    3.0941
    3.0827
    3.0690
    3.0685
    3.0701
    3.0618
    3.0502
    3.0297
    3.0245
    3.0091
    3.0039
    3.0062
    3.0058
    2.9793
    2.9620
    2.9689
    2.9489
    2.9459
    2.9476
    2.9301
    2.9112
    2.8999
    2.8960
    2.8883
    2.8899
    2.8675
    2.8693
    2.8512
    2.8478
    2.8400
    2.8432
    2.8189
    2.8083
    2.7990
    2.7788
    2.7888];
angle2FitVector = [-12.0006
    -11.7492
    -11.4996
    -11.2500
    -11.0004
    -10.7508
    -10.4994
    -10.2498
    -10.0002
    -9.7506
    -9.4992
    -9.2496
    -9.0000
    -8.7504
    -8.5008
    -8.2494
    -7.9998
    -7.7502
    -7.5006
    -7.2492
    -6.9996
    -6.7500
    -6.5004
    -6.2508
    -5.9994
    -5.7498
    -5.5002
    -5.2506
    -4.9992
    -4.7496
    -4.5000
    -4.2504
    -4.0008
    -3.7494
    -3.4998
    -3.2502
    -3.0006
    -2.7492
    -2.4996
    -2.2500
    -2.0004
    -1.7508
    -1.4994
    -1.2498
    -1.0002
    -0.7506
    -0.4992
    -0.2496
    0.0000
    0.2496
    0.4992
    0.7506
    1.0002
    1.2498
    1.4994
    1.7508
    2.0004
    2.2500
    2.4996
    2.7492
    3.0006
    3.2502
    3.4998
    3.7494
    4.0008
    4.2504
    4.5000
    4.7496
    4.9992
    5.2506
    5.5002
    5.7498
    5.9994
    6.2508
    6.5004
    6.7500
    6.9996
    7.2492
    7.5006
    7.7502
    7.9998
    8.2494
    8.5008
    8.7504
    9.0000
    9.2496
    9.4992
    9.7506
    10.0002
    10.2498
    10.4994
    10.7508
    11.0004
    11.2500
    11.4996
    11.7492
    12.0006
    12.0006
    11.7492
    11.4996
    11.2500
    11.0004
    10.7508
    10.4994
    10.2498
    10.0002
    9.7506
    9.4992
    9.2496
    9.0000
    8.7504
    8.5008
    8.2494
    7.9998
    7.7502
    7.5006
    7.2492
    6.9996
    6.7500
    6.5004
    6.2508
    5.9994
    5.7498
    5.5002
    5.2506
    4.9992
    4.7496
    4.5000
    4.2504
    4.0008
    3.7494
    3.4998
    3.2502
    3.0006
    2.7492
    2.4996
    2.2500
    2.0004
    1.7508
    1.4994
    1.2498
    1.0002
    0.7506
    0.4992
    0.2496
    0.0000
    -0.2496
    -0.4992
    -0.7506
    -1.0002
    -1.2498
    -1.4994
    -1.7508
    -2.0004
    -2.2500
    -2.4996
    -2.7492
    -3.0006
    -3.2502
    -3.4998
    -3.7494
    -4.0008
    -4.2504
    -4.5000
    -4.7496
    -4.9992
    -5.2506
    -5.5002
    -5.7498
    -5.9994
    -6.2508
    -6.5004
    -6.7500
    -6.9996
    -7.2492
    -7.5006
    -7.7502
    -7.9998
    -8.2494
    -8.5008
    -8.7504
    -9.0000
    -9.2496
    -9.4992
    -9.7506
    -10.0002
    -10.2498
    -10.4994
    -10.7508
    -11.0004
    -11.2500
    -11.4996
    -11.7492
    -12.0006];
% Perform linear fitting
p = polyfit(voltage2FitVector,angle2FitVector,1);
% Evaluate fitting with input voltage vector
angleVector = polyval(p,voltageVector);
end

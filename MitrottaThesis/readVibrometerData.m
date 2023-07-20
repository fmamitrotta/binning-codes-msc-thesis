function [timeVector,displacementVector] = readVibrometerData(filePath)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
% Open file
fileId = fopen(filePath);
% Acquire data
formatSpec = '%f %f';
dataArray = textscan(fileId,formatSpec,'HeaderLines',7,'Delimiter','\t');
timeVector = dataArray{1};
displacementVector = dataArray{2};
displacementVector = displacementVector-displacementVector(1);
% Close file
fclose(fileId);
end

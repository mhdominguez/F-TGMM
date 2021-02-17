function [xTrain, yTrain] = readTrainingDataBinary(filename)

fin = fopen(filename,'rb');
numFeatures = double( fread(fin,1,'int32') );
numSamples = double( fread(fin,1,'int32') );
xTrain = reshape( fread(fin, numFeatures * numSamples, 'float32'), [numFeatures numSamples] )';
yTrain = fread(fin, numSamples, 'float32');

fclose(fin);
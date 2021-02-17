function writeTrainingDataBinary(xTrain, yTrain, filename)

fin = fopen(filename,'wb');
numFeatures = size(xTrain,2);
numSamples = size(xTrain,1);
fwrite(fin,int32(numFeatures),'int32');
fwrite(fin,int32(numSamples),'int32');
xTrain = single(xTrain');
fwrite(fin, xTrain(:), 'float32') ;
fwrite(fin, single(yTrain(:)), 'float32');

fclose(fin);
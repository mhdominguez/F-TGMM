function generateCellDivisionDiscriminatorWithTemporalWindowClassifier(binaryFeatureFiles, numWeakLearners,numLeafs ,learnRate,classifierFileOut)



%numWeakLearners = 70;%30
%numLeafs = 60;%30
%learnRate = 0.1;

%classifierFileOut = 'C:\Users\Fernando\cppProjects\TrackingGaussianMixtures\NM2013-paperRun\src\filesToCopyToExeFolder\classifierCellDivisionWithTemporalWindow.txt';



%%
%read training data
xTrain = [];
yTrain = [];
for ii = 1:length(binaryFeatureFiles)
  [xTrainAux, yTrainAux] = readTrainingDataBinary(binaryFeatureFiles{ii});  
  xTrain = [ xTrain; xTrainAux];
  yTrain = [ yTrain; yTrainAux];
end

clear xTrainAux yTrainAux
   
disp(['Num. positive samples = ' num2str(sum(yTrain>=0.5)) '; Num. negative samples = ' num2str(sum(yTrain<0.5))]);
%%
%train classifier
rusTree = trainClassifierRusboost( xTrain, yTrain, numWeakLearners, numLeafs, learnRate);

%%
%save classifier
parseMatlabFitEnsembleToCppGentleBoostFormat(rusTree, classifierFileOut);


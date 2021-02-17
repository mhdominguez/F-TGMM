%run this code to generate a new classifier for the code
%function trainCellDivisionDiscriminatorWithTemporalWindow()

%{
%configuration for July 2014 classifier combining all datasets and annotations
numWeakLearners = 70;%30
numLeafs = 60;%30
learnRate = 0.1;
fileOut = 'C:\Users\Fernando\cppProjects\TrackingGaussianMixtures\NM2013-paperRun\src\filesToCopyToExeFolder\classifierCellDivisionWithTemporalWindow.txt';


trainingDatasets = {'E:\TGMMruns\GMEMtracking3D_2014_4_29_2_27_1_dataset12_08_28_drosophila_simview_temporalLRdeactivatedForCellDivisionTraining\annForCellDivDiscrWithTempWin\trainFeaturesCellDivDiscr_finalTraining_simmetry8.bin',...
                    'E:\TGMMruns\GMEMtracking3D_2014_4_29_2_27_1_dataset12_08_28_drosophila_simview_temporalLRdeactivatedForCellDivisionTraining\annForCellDivDiscrWithTempWin\trainFeaturesCellDivDiscr_finalTraining_symmetry8_subsampT2.bin',...
                    'E:\TGMMruns\GMEMtracking3D_2014_5_7_20_54_41_dataset12_09_24_zebrafish_simview_temporalLRdeactivatedForCellDivisionTraining\annForCellDivDiscrWithTempWin\trainFeaturesCellDivDiscr_finalTraining_symmetry8.bin',...
                    'E:\TGMMruns\GMEMtracking3D_2014_5_7_20_54_41_dataset12_09_24_zebrafish_simview_temporalLRdeactivatedForCellDivisionTraining\annForCellDivDiscrWithTempWin\trainFeaturesCellDivDiscr_finalTraining_symmetry8_subsampT2.bin',...
                    'E:\TGMMruns\GMEMtracking3D_2014_5_13_18_28_38_dataset12_08_28_drosophila_simview_CDTW_thr040\annForCellDivDiscrWithTempWin\trainFeaturesCellDivDiscr_finalTraining_symmetry8.bin',...
                    'E:\TGMMruns\GMEMtracking3D_2014_5_13_18_28_38_dataset12_08_28_drosophila_simview_CDTW_thr040\annForCellDivDiscrWithTempWin\trainFeaturesCellDivDiscr_finalTraining_symmetry8_subsampT2.bin',...
                    'E:\TGMMruns\GMEMtracking3D_2014_5_16_21_0_21_dataset12_10_09_zebrafish_confocal_CDTW_thr000_iter2\annForCellDivDiscrWithTempWin\trainFeaturesCellDivDiscr_trainingSet_symmetry8.bin'};
%}

%configuration for debugging new GUI in March 2015
numWeakLearners = 70;%30
minLeaf = 30;%30
learnRate = 0.1;

fileOut = './temp/classifierCDWT_MatlabScript.txt';

symmetry = 8;
trainingDatasets = {'Y:\Exchange\Fernando\NM2013-paperRun\TGMMruns\GMEMtracking3D_2015_1_27_17_45_20_Drosophila_trainCDWT_iter4\annForCellDivDiscrWithTempWin\trainFeaturesCellDivDiscr.bin'};

addpath c:\users\Fernando\matlabExternalToolboxes\boostingDemo-Torralba\matlabBoost
%%
%read training data
xTrain = [];
yTrain = [];
for ii = 1:length(trainingDatasets)
  [xTrainAux, yTrainAux] = readTrainingDataBinary(trainingDatasets{ii});  
  xTrain = [ xTrain; xTrainAux];
  yTrain = [ yTrain; yTrainAux];
end

clear xTrainAux yTrainAux
   
disp(['Num. positive samples = ' num2str(sum(yTrain>=0.5)) '; Num. negative samples = ' num2str(sum(yTrain<0.5))]);
%%
%train classifier
%rusTree = trainClassifierRusboost( xTrain, yTrain, numWeakLearners, numLeafs, learnRate);
figure;
[rusTree,thrVec,prec,rec] = trainClassifierRusboost( xTrain, yTrain, symmetry, numWeakLearners, minLeaf, learnRate, 5, gca);

%%
%save classifier
parseMatlabFitEnsembleToCppGentleBoostFormat(rusTree, fileOut);


%%
rmpath c:\users\Fernando\matlabExternalToolboxes\boostingDemo-Torralba\matlabBoost
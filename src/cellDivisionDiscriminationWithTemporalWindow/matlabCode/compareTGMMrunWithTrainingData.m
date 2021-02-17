function [FxTest, yTrain, tTrain] = compareTGMMrunWithTrainingData()

basenameTGMM = 'E:\TGMMruns\GMEMtracking3D_2015_3_6_17_42_12_Drosophila_12_08_28_trainCDWT_iter5\XML_finalResult_lht\GMEMfinalResult_frame';
annotationsFolder = 'Y:\Exchange\Fernando\NM2013-paperRun\TGMMruns\GMEMtracking3D_2015_1_27_17_45_20_Drosophila_trainCDWT_iter4\annForCellDivDiscrWithTempWin'
numTMdigits = 4;
anisotropyZ = 5.0;

%%
%read annotations
[objA, yTrain] = parseXmlAnnotationsFolder(annotationsFolder);

%extract time point for each annotation
tTrain = yTrain;
for ii = 1:length(yTrain)
   pp = findstr('TM', objA(ii).imFilename);
   pp = pp(1);
   tTrain(ii) =  str2double(objA(ii).imFilename(pp+2:pp+2+numTMdigits-1));
end


%%
%compare results for each time point
TM = unique(tTrain);
FxTest = yTrain; %contains the score in TGMM
for tt = 1:length(TM)
   pp = find(tTrain == TM(tt)); 
   Nt = length(pp);
   
   xyz = zeros(Nt, 3);
   for kk = 1:Nt
        xyz(kk,:) = objA(pp(kk)).m  .* [1 1 anisotropyZ];        
   end
   
   %read TGMM results
   addpath './cellDivAnnotatorGUI/readTGMM_XMLoutput/readTGMM_XMLoutput'   
   obj = readXMLmixtureGaussians([basenameTGMM num2str(TM(tt), '%.4d') '.xml']);
   rmpath './cellDivAnnotatorGUI/readTGMM_XMLoutput/readTGMM_XMLoutput'
   
   %parse xyz and cell division classifier score
   xyzTGMM = zeros(length(obj),3);
   Fx = zeros(length(obj),1);
   for kk = 1:length(obj)
       xyzTGMM(kk,:) = obj(kk).m .* [1 1 anisotropyZ];
       Fx(kk) = obj(kk).sigmaDistPrior;
   end
   
   %%
   %compare results
   [idx, dist] = knnsearch(xyzTGMM, xyz,'k',1);
   
   FxTest(pp) = Fx(idx);
   
end

%%
%display results
figure;hist(FxTest(yTrain < 0.5),[-1:0.05:1]);xlim([-1 1]);title('Classifier score for negative training samples');
figure;hist(FxTest(yTrain > 0.5),[-1:0.05:1]);xlim([-1 1]);title('Classifier score for positive training samples');



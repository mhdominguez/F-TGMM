function [FxTGMM, FxCpp, FxM, xTrainTGMM, xTrainCpp] = debugCheckFxScore(TM)

pathOutputTGMM = 'E:\TGMMruns\GMEMtracking3D_2015_3_13_19_38_11'
%TM = 35;
classifierFile = 'Y:\Exchange\Fernando\NM2013-paperRun\TGMMruns\GMEMtracking3D_2014_12_30_18_57_41_Drosophila_12_08_28_trainCDWT_noThr\classifierCellDivisionWithTemporalWindow_iter5_nW300_winRad5';

%%
%extract imgFilePattern and anisotropy
fid = fopen([pathOutputTGMM filesep 'experimentLog_0000.txt'],'r');
while 1
    tline = fgetl(fid);
    if ~ischar(tline), break, end
    
    if( length(tline) > 14 && strcmp(tline(1:14),'imgFilePattern') == 1)
        pos = strfind(tline,'=');
        imgFilePattern = tline(pos(1)+1:end);
    end
    
    if( length(tline) > 11 && strcmp(tline(1:11),'anisotropyZ') == 1)
        pos = strfind(tline,'=');
        anisotropyZ = str2double( tline(pos(1)+1:end) );        
    end
    
    if( length(tline) > 29 && strcmp(tline(1:29),'temporalWindowForLogicalRules') == 1)
        pos = strfind(tline,'=');
        temporalWindowForLogicalRules = tline(pos(1)+1:end);
    end      
end
fclose(fid);

%%
%read annotation
basenameTGMM = [pathOutputTGMM filesep 'XML_finalResult_lht' filesep 'GMEMfinalResult_frame'];
%read TGMM results
addpath './readTGMM_XMLoutput/readTGMM_XMLoutput'
obj = readXMLmixtureGaussians([basenameTGMM num2str(TM, '%.4d') '.xml']);
objCh = readXMLmixtureGaussians([basenameTGMM num2str(TM+1, '%.4d') '.xml']);
rmpath './readTGMM_XMLoutput/readTGMM_XMLoutput'


%parse children from objCh
ch = -ones(length(obj),2);
for ii = 1:length(objCh)
    parIdx = objCh(ii).parent + 1;
    if( parIdx > 0)
        if(ch(parIdx,1) < 0 )
            ch(parIdx,1) = ii; 
        else
            ch(parIdx,2) = ii; 
        end
    end
end

%find all the accepted cell divisions
xyz = zeros(length(obj),3);
FxTGMM = zeros(length(obj),1);
xyzN = 0;
for ii = 1:length(obj)
    if( ch(ii,2) < 0 )
        continue;
    end
    
    xyzN = xyzN + 1;
    xyz(xyzN,:) = obj(ii).m;
    FxTGMM(xyzN) = obj(ii).sigmaDistPrior;
end

xyz(xyzN+1:end,:) = [];
FxTGMM(xyzN+1:end) = [];

%%
%generate an annotation file to compute features
filenameAnnXML = [pathOutputTGMM filesep 'annForCellDivDiscrWithTempWin'];
mkdir(filenameAnnXML);
filenameAnnXML = [filenameAnnXML filesep 'debugAnnotations_TM' num2str(TM,'%.4d') '.xml'];


imgFilename = recoverFilenameFromPattern(imgFilePattern, TM);

fid = fopen(filenameAnnXML,'w');
fprintf(fid,'<?xml version="1.0" encoding="UTF-8"?>\n');
fprintf(fid,'<document>\n');

for ii = 1:size(xyz,1)
    %I only need the centroid to generate features (not even the file pattern)
    fprintf(fid,'<Surface name="Ellipsoid" id="1" numCoeffs="9" coeffs="-1 -1 -1 -1 -1 -1 %f %f %f" intensity="%f" covarianceMatrixSize="3" svFilename="empty" svIdx="1" imFilename="%s" class="cellDivisionWrong"></Surface>\n',xyz(ii,1),xyz(ii,2),xyz(ii,3), FxTGMM(ii), imgFilename);
end

fprintf(fid,'</document>\n');
fclose(fid);
 
 %%
 %run C++ code to generate features
%${TGMMrunRoot} <imgBasename> 5 symetry
TGMMrunRoot = regexprep(pathOutputTGMM,'\','/');

pathProg = [fileparts( mfilename('fullpath') ) filesep '..' filesep '..' filesep '..' filesep '..' filesep 'bin' filesep 'trainCellDivisionWithTemporalWindow'];
if( ispc )
    pathProg = [pathProg '.exe'];
end

cmd = [pathProg ' ' TGMMrunRoot ' ' imgFilePattern ' ' num2str(temporalWindowForLogicalRules) ' 1 1'];
[rr, ss] = system(cmd);

ss
if( rr ~= 0 )
    return;
end
 
 %%
 %load binary features
 binaryFeatureFiles = [pathOutputTGMM filesep 'annForCellDivDiscrWithTempWin' filesep 'trainFeaturesCellDivDiscr.bin'];
 [xTrain, yTrain] = readTrainingDataBinary(binaryFeatureFiles);  
 
 %load classifier
 load([classifierFile '.mat'], 'rusTree');
 %run classifier
 [label,Score] = predict(rusTree,xTrain);
 FxM = Score(:,2);
 
 %%
 %run C++ classifier using generated features and txt classifier
 classifierTxtFile = [classifierFile '.txt'];
 
 pathProg = [fileparts( mfilename('fullpath') ) filesep '..' filesep '..' filesep '..' filesep '..' filesep 'bin' filesep 'GentleBoost_PrecRecallCurve'];
 if( ispc )
     pathProg = [pathProg '.exe'];
 end
 
 cmd = [pathProg ' ' classifierTxtFile ' ' binaryFeatureFiles];
 [rr, ss] = system(cmd);
 
 ss
 if( rr ~= 0 )
     error 'Executring C++ code for precision recall curve'
 end
 
 %read results
 FxCpp = load([binaryFeatureFiles '_FXtraining.txt']);
 
%%
%load binary features from TGMM run-time
%you need to uncomment DEBUG_CDWT_FEATURES in TGMMsupport.cpp C++ code

if( nargout > 3 )
    binaryFeatureTGMM = [pathOutputTGMM filesep 'XML_finalResult_lht' filesep 'CDWTfeatures_TM' num2str(TM,'%.5d') '.bin'];
    xyzFeatureTGMM = [pathOutputTGMM filesep 'XML_finalResult_lht' filesep 'CDWTfeaturesXYZ_TM' num2str(TM,'%.5d') '.txt'];
    
    xyzTGMM = load(xyzFeatureTGMM);
    xTrainTGMM = readTrainingDataBinary(binaryFeatureTGMM);
    
    %find matching
    [idx, dist] = knnsearch(xyz,xyzTGMM, 'k', 1);
    
    idx = idx( dist < 1.0 );%everything might not have a match
    xTrainTGMM = xTrainTGMM( dist < 1.0, :);
    xTrainCpp = xTrain(idx,:);%so now they match
end
 
 %%
 %compare classifier scores
figure;
scatter(FxM,FxTGMM)
xlabel('Matlab Fx score')
ylabel('TGMM C++ Fx score')
grid on
title(['Frame ' num2str(TM)])

figure;
scatter(FxCpp,FxTGMM)
xlabel('C++ Fx score')
ylabel('TGMM C++ Fx score')
grid on
title(['Frame ' num2str(TM)])


figure;
scatter(FxCpp,FxM)
xlabel('C++ Fx score')
ylabel('Matlab Fx score')
grid on
title(['Frame ' num2str(TM)])
 
 
 
 
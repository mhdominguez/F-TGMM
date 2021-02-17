%script to resave Drosophila annotation in order to use them with new CDTW interface

dataset = 2;%1->12-08-28; 2->12-09-24 (zebrfish)

%%
switch(dataset)
    case 1
        %list of folders containing Drosophila cell division annotations
        folderAnn = {'Y:\Exchange\Fernando\NM2013-paperRun\TGMMruns\GMEMtracking3D_2014_12_30_18_57_41_Drosophila_12_08_28_trainCDWT_noThr',...
            'Y:\Exchange\Fernando\NM2013-paperRun\TGMMruns\GMEMtracking3D_2015_1_2_4_50_40_Drosophila_12_08_28_trainCDWT_iter2',...
            'Y:\Exchange\Fernando\NM2013-paperRun\TGMMruns\GMEMtracking3D_2015_1_27_17_45_20_Drosophila_trainCDWT_iter4',...
            'E:\TGMMruns\GMEMtracking3D_2015_3_6_17_42_12_Drosophila_12_08_28_trainCDWT_iter5',...
            };
        
        %'Y:\Exchange\Fernando\NM2013-paperRun\TGMMruns\GMEMtracking3D_2015_1_2_18_35_49_Drosophila_12_08_28_trainCDWT_iter3',... contains annotations for test TM
        TMpos = [154:157];
        %TGMM folder run twith the new CDTW code
        outputTGMMfolder = 'E:\TGMMruns\GMEMtracking3D_2015_3_16_18_21_34_Drosophila_12_08_28_trainCDWT_V2_iter0'
        
    case 2
        %list of folders containing Drosophila cell division annotations
        folderAnn = {'Y:\Exchange\Fernando\NM2013-paperRun\TGMMruns\GMEMtracking3D_2014_12_31_3_33_3_Zebrafish_12_09_24_trainCDWT_noThr',...
            'Y:\Exchange\Fernando\NM2013-paperRun\TGMMruns\GMEMtracking3D_2015_1_10_4_10_13_Zebrafish_12_09_24_trainCDWT_iter2',...            
            };
        
        %'Y:\Exchange\Fernando\NM2013-paperRun\TGMMruns\GMEMtracking3D_2015_1_2_18_35_49_Drosophila_12_08_28_trainCDWT_iter3',... contains annotations for test TM
        TMpos = [119:122];
        %TGMM folder run twith the new CDTW code
        outputTGMMfolder = 'E:\TGMMruns\GMEMtracking3D_2015_3_17_20_30_20_Zebrafish_12_09_24_trainCDWT_v2_iter0'
        
end
        
%%
%read annotations
obj = [];
for ii = 1:length(folderAnn)
    [objAux, ff] = parseXmlAnnotationsFolder([folderAnn{ii} filesep 'annForCellDivDiscrWithTempWin']);
    obj = [obj; objAux];
end

clear objAux;

%%
%parse xyzt for each annotation
annN = length(obj);
xyzt = zeros(annN, 4);

for ii = 1:annN
    xyzt(ii,:) = [obj(ii).m, str2double(obj(ii).imFilename(TMpos))];
end

%%
%find possible match for each annotation
addpath 'cellDivAnnotatorGUI\readTGMM_XMLoutput\readTGMM_XMLoutput'
erase = [];
for TM = unique(xyzt(:,4))'
        
    pos = find(xyzt(:,4) == TM);
    
    %readXML folder
    objTGMM = readXMLmixtureGaussians([outputTGMMfolder filesep 'XML_finalResult_lht' filesep 'GMEMfinalResult_frame' num2str(TM,'%.4d') '.xml']);
    %keep only cell divisions
    objTGMMCh = readXMLmixtureGaussians([outputTGMMfolder filesep 'XML_finalResult_lht' filesep 'GMEMfinalResult_frame' num2str(TM + 1,'%.4d') '.xml']);
    
    ch = -ones(length(objTGMM), 2);
    for jj = 1:length(objTGMMCh)
        parIdx = objTGMMCh(jj).parent + 1;
        if( parIdx > 0 )
           if( ch(parIdx,1) < 0 )
               ch(parIdx,1) = jj;
           else
               ch(parIdx,2) = jj;
           end
        end
    end
    
    %read centroid
    xyzTGMM = zeros(length(objTGMM), 3);    
    for jj = 1:size(xyzTGMM,1)
        xyzTGMM(jj,:) = objTGMM(jj).m;
    end
      
    
    %find match
    [idx, dist] = knnsearch(xyzTGMM, xyzt(pos,1:3), 'k',1);
    erase = [erase; pos( dist > 1.0 | ch(idx,2) <= 0 )];
end
rmpath 'cellDivAnnotatorGUI\readTGMM_XMLoutput\readTGMM_XMLoutput'

%%
%save matches
disp(['Matched ' num2str(length(obj)-length(erase)) ' out of ' num2str(length(obj)) ' annotations']);
obj(erase) = [];
addpath 'cellDivAnnotatorGUI'

%save all the annotations as an XML surface object
%no header or footer to be able to combine multiple savings with CAT
%easily
pathAnn = [outputTGMMfolder filesep 'annForCellDivDiscrWithTempWin' filesep];
fileOut = [pathAnn 'classifierAnnotations_' datestr(now,30) '.xml'];
fout=fopen(fileOut,'w');
%write header
fprintf(fout,'<?xml version="1.0" encoding="UTF-8"?>\n<document>\n');
%write body
numCoeffs = 9;
for kk = 1:length(obj)
    blob = obj(kk);
    fprintf(fout,'<Surface name="Ellipsoid" id="1" numCoeffs="%d" ',numCoeffs);
    fprintf(fout,'coeffs="');
    coeffs = [blob.W([1 2 3 5 6 9]) blob.m];
    fprintf(fout,'%f ',coeffs);
    
    
    fprintf(fout,'" intensity="%g',-1.0);
    
    
    fprintf(fout,'" covarianceMatrixSize="%d" \n',3);
    fprintf(fout,'imFilename="%s" class="%s"></Surface>\n',blob.imFilename,blob.class);    
end
%write footer
fprintf(fout,'</document>\n');
fclose(fout);

rmpath 'cellDivAnnotatorGUI'

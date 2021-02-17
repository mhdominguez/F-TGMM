function [thrVec,prec,rec] = fitMidplaneThr(midPlaneVal, cellDivLabel,Kcv)

posP = find(cellDivLabel > 0.5);
posPN = length(posP);
posN = find(cellDivLabel <= 0.5);
posNN = length(posN);



kcvP = crossvalind('Kfold', posPN, Kcv); 
kcvN = crossvalind('Kfold', posNN, Kcv); 

%calculate range of thresholds
thrMin = min(midPlaneVal(posN)) - 0.1;
thrMax = max(midPlaneVal(posP)) + 0.1;
thrVec = thrMin:0.05:thrMax;

prec = zeros(length(thrVec),Kcv);
rec = zeros(length(thrVec),Kcv);
for ii = 1:Kcv

    xTest = [midPlaneVal(posP(kcvP==ii)); midPlaneVal(posN(kcvN==ii))];
    yTest = [cellDivLabel(posP(kcvP==ii)); cellDivLabel(posN(kcvN==ii))];
    
    for kk = 1:length(thrVec)    
        thr = thrVec(kk);
        
        yPred = xTest < thr;
        
        TP = sum(yPred > 0.5 & yTest > 0.5);
        TN = sum(yPred <= 0.5 & yTest <= 0.5);
        FP = sum(yPred > 0.5 & yTest <= 0.5);
        FN = sum(yPred <= 0.5 & yTest > 0.5);
        
        prec(kk,ii) = TP / (TP+FP);
        rec(kk,ii) = TP / (TP+FN);
    end
end

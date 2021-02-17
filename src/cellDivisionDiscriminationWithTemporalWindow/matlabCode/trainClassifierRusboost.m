%following example from http://www.mathworks.com/help/stats/ensemble-methods.html#btgw1m1
%for imbalanced data

%we do not separate training data in test set and train set. This is for
%final trainings. No cross-validation either

%use parseMatlabFitEnsembleToCppGentleBoostFormat(rusTree, filenameOut) to save rusTree model
%as a text file that can be used in our C/C++ code
function  [rusTree,thrVec,prec,rec] = trainClassifierRusboost( xTrainOrig, yTrainOrig, symmetry, numWeakLearners, minLeaf, learnRate, Kcv, axes2)

displayTrainingError = false;%set to true for debugging purposes



disp 'Tabulate all data before cross validation'
tabulate( yTrainOrig )

disp 'Learning classifier'
t = ClassificationTree.template('minleaf',minLeaf);%sum(yTrain(istrain)>0) for stumps if RatioToSmallest is [1 1]; otherwise 5 is a good value for large trees

if( displayTrainingError )
    disp '=====WARNING:displaying training error=============='
end
invMap = [];%parfor needs it
xTrainFull = [];
yTrainFull = [];
%-----------------------------------------
if( symmetry > 1 )
    %account for symmetry repetition in the training set (so it would be too
    %optimistic)
    
    nReal = length(yTrainOrig) / symmetry; 
    
    fRepeat = -1;
    for ii = 1:size(xTrainOrig,2)
        if( nReal == length( unique(xTrainOrig(:,ii)) ) )
            fRepeat = ii;%use this column (feature) to identify repetition
           break; 
        end
    end
    if( fRepeat == -1)%sometimes we have repeats and we have to use th emode
        qq = zeros(size(xTrainOrig,2),1);
        for ii = 1:size(xTrainOrig,2)
            qq(ii) = length( unique(xTrainOrig(:,ii)) );            
        end
        fRepeat = find(qq == mode(qq));
        fRepeat = fRepeat(1);
    end
    
    
    
    if( fRepeat > 0 )       
        %length(IA) = nReal; xTrainOrig(IA,:) contains unique samples
        [~,IA,IC] = unique(xTrainOrig(:,fRepeat));%length(IC) = size(xTrainFull,2) and find(IC == ii) = symmetry for ii = 1:nReal        
    end
    
    %make sure they all have symmetry copies
    [u, p] = hist(IC,[min(IC):max(IC)]);    
    ppp = find(u ~= symmetry);
    erase = [];
    for ii = 1:length(ppp)
       erase = [erase; find(IC == ppp(ii))];
    end
    
    if( isempty(erase) == false )
       yTrainOrig(erase) = [];
       xTrainOrig(erase,:) = [];
       
       [~,IA,IC] = unique(xTrainOrig(:,fRepeat));%length(IC) = size(xTrainFull,2) and find(IC == ii) = symmetry for ii = 1:nReal   
        nReal = length(yTrainOrig) / symmetry;
    end
    
    
    
    %copy before restoring original samples
    xTrainFull = xTrainOrig;
    yTrainFull = yTrainOrig;
    
    xTrainOrig = xTrainOrig(IA,:);
    yTrainOrig = yTrainOrig(IA);
    
    
    %generate "invers map" from IC (kind of a multimap container)
    invMap = zeros(nReal,symmetry);
    invMapN = zeros(nReal,1);%counter
    for ii = 1:length(IC)
       row = IC(ii);
       col = invMapN(row) + 1;
       invMapN(row) = invMapN(row) + 1;
       invMap(row, col) = ii;
    end

end


%------------------------------------------

posP = find(yTrainOrig > 0.5);
posPN = length(posP);
posN = find(yTrainOrig <= 0.5);
posNN = length(posN);
numFeatures = size(xTrainOrig,2);


kcvP = crossvalind('Kfold', posPN, Kcv); 
kcvN = crossvalind('Kfold', posNN, Kcv); 


rusTree = cell(Kcv,1);
ppCell = cell(Kcv,1);
qqCell = cell(Kcv,1);
plotArray = zeros(Kcv,numWeakLearners);

if( displayTrainingError )
    yKcv = NaN(size(yTrainFull));
else
    yKcv = NaN(size(yTrainOrig));
end

if( matlabpool('size') > 0 && matlabpool('size') ~= Kcv )
    matlabpool('close');
end
if( matlabpool('size') == 0 )
    matlabpool(Kcv);
end

%disp '=======WARNING: PARFOR DISABLED======'
parfor ii = 1:Kcv

    xTest = [xTrainOrig(posP(kcvP==ii),:); xTrainOrig(posN(kcvN==ii),:)];
    yTest = [yTrainOrig(posP(kcvP==ii)); yTrainOrig(posN(kcvN==ii))];
    
    %generate extended features
    if( symmetry > 1)
        A = ([posP(kcvP~=ii);posN(kcvN~=ii)]);
        
        pp = invMap(A,:);
        pp = pp(:);
        xTrain = xTrainFull(pp,:);
        yTrain = yTrainFull(pp);
        
    else%we can use all the features directly
        xTrain = [xTrainOrig(posP(kcvP~=ii),:); xTrainOrig(posN(kcvN~=ii),:)];
        yTrain = [yTrainOrig(posP(kcvP~=ii)); yTrainOrig(posN(kcvN~=ii))];
    end
    rusTree{ii} = fitensemble(xTrain,yTrain,'RUSBoost',numWeakLearners,t,'LearnRate',learnRate);
    
    
    if( displayTrainingError )
        xTest = xTrain;
        yTest = yTrain;
    end
    
    %obtain scores from C++ program (so it exactly matches what we do in
    %the tracking code)
    
    basename = [tempname '_' num2str(ii)];
    filenameClassifier = [basename '.txt'];
    filenameFeatures = [basename '.bin'];
    
    %write features
    fin = fopen(filenameFeatures,'wb');
    fwrite(fin,int32(numFeatures),'int32');
    fwrite(fin,length(yTest),'int32');
    
    xTest = xTest';
    fwrite(fin,single(xTest(:)),'float32');
    fwrite(fin,single(yTest),'float32');        
    fclose(fin);
    xTest = xTest';
    
    %write classifier
    parseMatlabFitEnsembleToCppGentleBoostFormat(rusTree{ii}, filenameClassifier);
    
    
    %run C++ code
    pathProg = [fileparts( mfilename('fullpath') ) filesep '..' filesep '..' filesep '..' filesep  'bin' filesep 'GentleBoost_PrecRecallCurve'];
    if( ispc )
        pathProg = [pathProg '.exe'];
    end
    
    cmd = [pathProg ' ' filenameClassifier ' ' filenameFeatures];    
    [rr, ss] = system(cmd);
    
    if( rr ~= 0 )
        ss
        error 'Executring C++ code for precision recall curve'
    end
    
    %read results
    qqCell{ii} = load([filenameFeatures '_FXtraining.txt']);
    if( displayTrainingError )
        ppCell{ii} = pp;
    else
        ppCell{ii} = [posP(kcvP==ii) ; posN(kcvN==ii)];
    end
    %yKcv(pp) = qq;  
    
    
    plotArray(ii,:) = loss(rusTree{ii},xTest,yTest,'mode','cumulative')';
    
    
    
    
end


if( isempty(axes2) == false)   
    cla(axes2);
    hold(axes2,'on');
    for ii = 1:Kcv
        plot(axes2, plotArray');
    end
    title(axes2,['Using ' num2str(Kcv) '-fold cross validation']);
    xlabel(axes2,'Number of weak classifiers');
    ylabel(axes2,'Error rate');
    hold(axes2,'off');
    grid(axes2,'on');    
end

%reduce results
for ii = 1:Kcv
   yKcv(ppCell{ii}) = qqCell{ii}; 
end

%calculate precision recall
thrVec = 0:0.02:1.0;
prec = zeros(length(thrVec),1);
rec = zeros(length(thrVec),1);

if( displayTrainingError  && symmetry > 1)
    yTrainOrig = yTrainFull;
    xTrainOrig = xTrainFull;
end

for kk = 1:length(thrVec)
    thr = thrVec(kk);
    
    yPred = yKcv > thr;
    
    TP = sum(yPred > 0.5 & yTrainOrig > 0.5);
    TN = sum(yPred <= 0.5 & yTrainOrig <= 0.5);
    FP = sum(yPred > 0.5 & yTrainOrig <= 0.5);
    FN = sum(yPred <= 0.5 & yTrainOrig > 0.5);
    
    prec(kk) = TP / (TP+FP);
    rec(kk) = TP / (TP+FN);
end
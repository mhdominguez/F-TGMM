posP = find( yTestZC > 0.5 );
nf = 50;%number of features to calculate distance
knn = 10;


ww =  imp(1:nf) / sum(imp(1:nf));%normalized weights
[fidx, fdist]  = knnsearch( bsxfun( @times, xTestZS(:,impIdx(1:nf)), ww), bsxfun(@times, xTestZC(posP,impIdx(1:nf)), ww), 'k', knn );


%%
%
ii = 6;

[boxCell, label, frameVec] = debugLoadBoxesInCDTW(posP(ii)-1, 'E:\temp\3DHaarBoxes\zebrafish_confocal', true);
pathBoxes = 'E:\temp\3DHaarBoxes\drosophila_simview';
for jj = 1:size(fidx,2)
    [boxCell, label, frameVec] = debugLoadBoxesInCDTW(fidx(ii,jj)-1, pathBoxes, true);

end

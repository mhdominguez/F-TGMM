pathBoxes = 'E:\temp\3DHaarBoxes';
basename = 'drosophila_simview';

ii = 0;
[boxCell, label, frameVec] = debugLoadBoxesInCDTW(ii, [pathBoxes '\' basename], true);
while( isempty( boxCell ) == false )
    saveas(gcf,['./png/' basename '_y_' num2str(label) '_' num2str(ii,'%.4d') '.png'])
    close;
    ii = ii + 1;
    [boxCell, label, frameVec] = debugLoadBoxesInCDTW(ii, [pathBoxes '\' basename], true);
end
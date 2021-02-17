pathBoxes = 'E:\temp\3DHaarBoxes';
%basename = 'drosophila_simview';
basename = 'zebrafish_confocal';

ii = 0;
intensity = [];
boxCell = 1;
while( isempty( boxCell ) == false )        
    [boxCell, label, frameVec] = debugLoadBoxesInCDTW(ii, [pathBoxes '\' basename], false);
    
    if( label > 0.5 )
        qq = cell2mat(boxCell);
        intensity = [intensity; qq(:)];
    end
    
    ii = ii + 1;
end
%in coordination with c++ code debug_mainSingleWindowForDaughters_writeImageBoxes()
function [boxCell, label, frameVec] = debugLoadBoxesInCDTW(sampleNumber, pathBoxes, debugDisplay)

%pathBoxes = 'E:\temp\3DHaarBoxes';

fileTxt = [pathBoxes '\aaa_boxIndex.txt'];%contains the id of each binary box
basenameBox = [pathBoxes '\box_'];%basename to each box

%open text file with all the info
%columns are boxId (for bin file), sample id, time point, yTrain value
boxIdx = load(fileTxt);


%load all the boxes belonging to the sample number
pos = find(boxIdx(:,2) == sampleNumber);
N = length(pos);

if( isempty(pos) )
   boxCell = [];
   label = nan(1);
   frameVec = [];
   return;
end

boxCell = cell(N,1);

for ii = 1:N
    boxName = [basenameBox num2str(boxIdx(pos(ii),1),'%.6d') '.bin'];
    fid = fopen(boxName, 'rb');
    
    %read dimensions
    boxRadius = fread(fid,3,'int32');
    boxSize = 1+ 2 * boxRadius';
    boxCell{ii} = reshape( fread(fid,prod(boxSize),'float32'), boxSize );    
    fclose(fid);    
end

label = boxIdx(pos(ii),4);

%reorder based on time
[frameVec,idx] = sort(boxIdx(pos,3), 'ascend');
boxCell = boxCell(idx); 

%%
if( debugDisplay )
    %find normalizing coefficients
    qq = cell2mat(boxCell);
    thrI = prctile(qq(:),[3 97]);
    %display elements
    figure;
    for ii = 1:length(boxCell)
        box = boxCell{ii};
        
        subplot(3,4,ii);
        cc = (size(box,3) + 1)/2;
        %im = squeeze(box(:,:,cc));%midplane
        im = max(box(:,:,cc-1:cc+1),[],3);%MIP of a slab
        
        im = (im - thrI(1)) / (thrI(2)-thrI(1));
        im( im<0 ) = 0;
        im( im > 1 ) = 1;
        imshow(im);
        title(['TM = ' num2str(frameVec(ii)) ';y = ' num2str(label) ]);
    end
    
end

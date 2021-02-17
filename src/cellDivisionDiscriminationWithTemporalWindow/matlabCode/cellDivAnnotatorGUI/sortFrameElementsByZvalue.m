%sorts elements of a frame by size

function [sortIdxList, aux]=sortFrameElementsByZvalue(blobStructFrame)



vol=zeros(length(blobStructFrame),1);
for ii=1:length(blobStructFrame)
    
    if(isempty(blobStructFrame(ii).frame))
        vol(ii)=1e32;
    else
        vol(ii) = blobStructFrame(ii).center(3);        
    end
end

[aux sortIdxList]=sort(vol,'ascend');%returns on top the elements with highest possibility of being background


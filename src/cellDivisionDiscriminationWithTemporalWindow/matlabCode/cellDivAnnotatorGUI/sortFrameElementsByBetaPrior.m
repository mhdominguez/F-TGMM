%sorts elements of a frame by size

function [sortIdxList, aux]=sortFrameElementsByBetaPrior(blobStructFrame)



vol=zeros(length(blobStructFrame),1);
for ii=1:length(blobStructFrame)
    
    if(isempty(blobStructFrame(ii).frame))
        vol(ii)=-1e32;
    else
        vol(ii) = blobStructFrame(ii).radius;        
    end
end

[aux sortIdxList]=sort(vol,'descend');%returns on top the elements with highest possibility of being background


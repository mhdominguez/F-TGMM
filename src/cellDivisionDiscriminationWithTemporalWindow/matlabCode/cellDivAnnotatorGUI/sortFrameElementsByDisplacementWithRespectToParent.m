%sorts elements of a frame by size

function [sortIdxList, aux ]=sortFrameElementsByDisplacementWithRespectToParent(blobStructFrame,blobStructFrameChildren,anisotropy)


vol=zeros(length(blobStructFrame),1);
for ii=1:length(blobStructFrame)
    
    
    if(blobStructFrame(ii,3)<-1e31)%dead cell
        vol(ii)=0;
    else
        ch=blobStructFrame(ii,17);
        
        if(ch >= 0)
            if(blobStructFrameChildren(ch+1,3)>-1)%parent not dead
                vol(ii)=norm((blobStructFrameChildren(ch+1,3:5)-blobStructFrame(ii,3:5)).*[1 1 anisotropy]);
            end
        end
    end
    
end


[aux, sortIdxList]=sort(vol,'descend');%returns 
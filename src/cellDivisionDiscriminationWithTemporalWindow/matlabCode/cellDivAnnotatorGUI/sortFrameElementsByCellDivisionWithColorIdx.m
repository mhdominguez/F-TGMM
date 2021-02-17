%sort elements weather if there are cell divisions or not. We need
%.colorIdx on the structure

function sortIdxList=sortFrameElementsByCellDivisionWithColorIdx(blobStructFrame)



vol=zeros(length(blobStructFrame),1);
for ii=1:length(blobStructFrame)
    
    if(isempty(blobStructFrame(ii).frame))
        vol(ii)=0;
    else
        
        
        
        if(blobStructFrame(ii).surface.coeffs(7)<-1e30)%dead cell
            vol(ii)=0;
        elseif(blobStructFrame(ii).colorIdx==5)
            vol(ii)=1;
        end
    end
end

[aux sortIdxList]=sort(vol,'descend');%returns 

%randomized the order of teh cell division
pos=find(aux==1);
sortIdxList(pos)=sortIdxList(pos(randperm(length(pos))));
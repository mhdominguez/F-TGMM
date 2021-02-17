function [sortIdxList, aux]=sortFrameElementsByTrackBirth(blobStructFrame)

vol=zeros(length(blobStructFrame),1);
for ii=1:length(blobStructFrame)
    
    if(isempty(blobStructFrame(ii).frame))
        vol(ii)=-1e32;
    else
        
        
        
        if(blobStructFrame(ii).surface.coeffs(7)<-1e30)%dead cell
            vol(ii)=-1e32;
        elseif(length(blobStructFrame(ii).solutions(1).childrenIdx)>2)
            %make sure both children are alive
            par=blobStructFrame(ii).solutions(1).parentIdx;
            if( par(1) > 4.29e9 )
                vol(ii) = 0;
            else
                vol(ii) = -1;
            end
        end
    end
end

[aux sortIdxList]=sort(vol,'descend');%returns 

% % % % %randomized the order of teh cell division
% % % % pos=find(aux==1);
% % % % sortIdxList(pos)=sortIdxList(pos(randperm(length(pos))));
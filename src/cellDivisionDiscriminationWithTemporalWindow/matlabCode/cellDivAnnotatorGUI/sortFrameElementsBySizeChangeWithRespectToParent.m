%sorts elements of a frame by size

function sortIdxList=sortFrameElementsBySizeChangeWithRespectToParent(blobStructFrame,blobStructFrameChildren,anisotropy)

scale = [1 1 anisotropy];
scale = scale' * scale;

vol=zeros(length(blobStructFrame),1);
sigma1=zeros(3);
for ii=1:length(blobStructFrame)
    
    if(isempty(blobStructFrame(ii).frame))
        vol(ii)=0;
    else
        if(blobStructFrame(ii).center(1)<-1e31)%dead cell
            vol(ii)=0;
        else
            ch=blobStructFrame(ii).solutions(1).childrenIdx;
            
            if(length(ch) == 2 )
                if(blobStructFrameChildren(ch(2)+1).center(1)>-1)%parent not dead
                    coeffs=blobStructFrame(ii).surface.coeffs;
                    sigma1(1,1)=coeffs(1);sigma1(2,1)=coeffs(2);sigma1(3,1)=coeffs(3);
                    sigma1(2,2)=coeffs(4);sigma1(3,2)=coeffs(5);sigma1(3,3)=coeffs(6);
                    sigma1(1,2)=sigma1(2,1);sigma1(1,3)=sigma1(3,1);sigma1(2,3)=sigma1(3,2);                    
                    WPar=sigma1./scale;
                    
                    
                    coeffs=blobStructFrameChildren(ch(2)+1).surface.coeffs;                    
                    sigma1(1,1)=coeffs(1);sigma1(2,1)=coeffs(2);sigma1(3,1)=coeffs(3);
                    sigma1(2,2)=coeffs(4);sigma1(3,2)=coeffs(5);sigma1(3,3)=coeffs(6);
                    sigma1(1,2)=sigma1(2,1);sigma1(1,3)=sigma1(3,1);sigma1(2,3)=sigma1(3,2);                    
                    WCh=sigma1./scale;
                    
                    vol(ii)= abs( det(WPar) - det(WCh) ) / det(WPar);
                end
            end
        end
    end
end


[aux sortIdxList]=sort(vol,'descend');%returns from larger changes in size
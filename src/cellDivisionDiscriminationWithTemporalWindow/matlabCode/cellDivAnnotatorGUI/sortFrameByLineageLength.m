%checks the length of the lineage

function [sortIdxList, sortValList]=sortFrameByLineageLength(frameId)

global blobStructGlobal;

vol=zeros(length(blobStructGlobal(frameId + 1,:)),1);

for ii=1:length(vol)
    
    if(isempty(blobStructGlobal(frameId+1,ii).frame))
        vol(ii)=1e32;%they will be displayed teh last
    else
        
        %navigate the tree upstream
        ll = 1;
        ff = frameId + 1;
        bb = ii;
        par = blobStructGlobal(ff,bb).solutions.parentIdx;        
        while( par(1) < 4294967290 )
            ll = ll +1;
            ff = par(1) + 1;
            bb = par(2) + 1;
            par = blobStructGlobal(ff,bb).solutions.parentIdx;
        end
        
        %navigate tree dowstream (we choose always the left child)
        ff = frameId + 1;
        bb = ii;
        ch = blobStructGlobal(ff,bb).solutions.childrenIdx;        
        while( isempty(ch) == false )
            ll = ll +1;
            ff = ch(1) + 1;
            bb = ch(2) + 1;
            ch = blobStructGlobal(ff,bb).solutions.childrenIdx;
        end
                
        
        vol(ii) = ll;
    end
end

[sortValList, sortIdxList]=sort(vol,'ascend');%we will see first the shortest lineages
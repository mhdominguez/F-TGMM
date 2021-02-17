%sorts elements of a frame by number of supervoxels asosciated with them

function [sortIdxList, aux]=sortFrameElementsByNumberOfSupervoxels(blobStructFrame)

%{
if(matlabpool('size')~=8)
    if(matlabpool('size')>0)
        matlabpool close;
    end
    matlabpool(8);
end
%}

vol=zeros(length(blobStructFrame),1);
for ii=1:length(blobStructFrame)
    
    if(isempty(blobStructFrame(ii).frame))
        vol(ii)=-1;
    else
        
        if(blobStructFrame(ii).center(1)<-1e30 )%dead cell
            vol(ii) = -1;
        else
            vol(ii) = length(blobStructFrame(ii).svIdx);
        end
    end
end

vol(vol<eps)=eps;

[aux, sortIdxList]=sort(vol,'descend');%returns 
%sorts elements of a frame by size

function [sortIdxList aux]=sortFrameElementsBySize(blobStructFrame,scale)

%{
if(matlabpool('size')~=8)
    if(matlabpool('size')>0)
        matlabpool close;
    end
    matlabpool(8);
end
%}

auxScale=scale'*scale;
vol=zeros(length(blobStructFrame),1);
for ii=1:length(blobStructFrame)
    
    if(isempty(blobStructFrame(ii).frame))
        vol(ii)=1e32;
    else
        coeffs=blobStructFrame(ii).surface.coeffs;
        sigma1=zeros(3);
        sigma1(1,1)=coeffs(1);sigma1(2,1)=coeffs(2);sigma1(3,1)=coeffs(3);
        sigma1(2,2)=coeffs(4);sigma1(3,2)=coeffs(5);sigma1(3,3)=coeffs(6);
        sigma1(1,2)=sigma1(2,1);sigma1(1,3)=sigma1(3,1);sigma1(2,3)=sigma1(3,2);
        
        sigma1=sigma1./auxScale;
        
        %volume of ellipse is 4/3*pi*a*b*c. a=1/sqrt(lamba_1) -> vol is
        %proportional to 1./sqrt(det(sigma1)) -> just calculated det(sigma1)
        %and sort in ascending order (equivalent to returning from larger size to smaller)
        if(coeffs(7)<-1e30 || blobStructFrame(ii).center(1)<-1e30)%dead cell
            vol(ii)=1e32;
        else
            vol(ii)=det(sigma1);
            %vol(ii)=min([det(sigma1([1 2],[1 2])) det(sigma1([2 3],[2 3])) det(sigma1([1 3],[1 3]))]) ;
        end
    end
end

vol(vol<eps)=eps;

[aux sortIdxList]=sort(vol,'ascend');%returns 
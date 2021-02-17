%sorts elements of a frame by size

function sortIdxList=sortFrameElementsByEccentricity(blobStructFrame,scale)

%{
if(matlabpool('size')~=8)
    if(matlabpool('size')>0)
        matlabpool close;
    end
    matlabpool(8);
end
%}

vol=zeros(length(blobStructFrame),1);

auxScale=scale'*scale;
for ii=1:length(blobStructFrame)
    
    if(isempty(blobStructFrame(ii).frame))
        vol(ii)=0;
    else
        coeffs=blobStructFrame(ii).surface.coeffs;
        sigma1=zeros(3);
        sigma1(1,1)=coeffs(1);sigma1(2,1)=coeffs(2);sigma1(3,1)=coeffs(3);
        sigma1(2,2)=coeffs(4);sigma1(3,2)=coeffs(5);sigma1(3,3)=coeffs(6);
        sigma1(1,2)=sigma1(2,1);sigma1(1,3)=sigma1(3,1);sigma1(2,3)=sigma1(3,2);
        
        sigma1=sigma1./auxScale;
        
        if(coeffs(7)<-1e30 || blobStructFrame(ii).center(1)<-1e30)%dead cell
            vol(ii)=0;
        else
            [V D]=eig(sigma1);
            D=diag(D);
            if(min(D)<1e-8)
                vol(ii)=1e32;
            else
                vol(ii)=max(D)/min(D);
            end
        end
    end
end

vol(vol<eps)=eps;

[aux sortIdxList]=sort(vol,'descend');%returns 
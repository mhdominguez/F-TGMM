%displays lineage for all neighbors around a given point

function [rectCoord Zstep]=neighboringLineageVisualization(frameOffset,blobStruct,frameId,blobId,numNN,deltaT,anisotropy,h,pyramidLevelValue)

Zstep=30;%space between time points to mltiplex the Z axis between space and time

frameAux=frameId-frameOffset;%to avoid passing all blobStruct since we are only looking locally

if(isempty(blobStruct(frameAux,blobId).frame))
    disp 'Selected blob is non-existant'
    return;
end

if(blobStruct(frameAux,blobId).center(1)<-1e31 || blobStruct(frameAux,blobId).surface.coeffs(7)<-1e31)
    disp 'Selected blob is dead'
    return;
end


%parse xyz components for given time point
N=size(blobStruct,2);
xyz=zeros(N,3);

for kk=1:N
   if(isempty(blobStruct(frameAux,kk).frame))
       break;
   end
   xyz(kk,:)=blobStruct(frameAux,kk).center;
end

xyz(kk:end,:)=[];
xyz(:,3)=xyz(:,3)*anisotropy;

%find nearest neighbors
xyzCenter=xyz(blobId,:);
dd=sum((xyz-repmat(xyzCenter,[size(xyz,1) 1])).^2,2);
[ddSort ddIdx]=sort(dd,'ascend');

kNN=ddIdx(1:numNN+1);%the first one is always itself

hold(h,'on');
%-----------------------------------------

rectCoord=[1e32 0 1e32 0];%xMin xMax yMin yMax
for kk=1:length(kNN)
   rectCoordAux=displayPartialLineageTree(Zstep,frameOffset,blobStruct,frameId-1,kNN(kk)-1,deltaT,1,anisotropy,h,kk,pyramidLevelValue);    
   disp(['frameId=' num2str(frameId-1) ';blobId=' num2str(kNN(kk)-1) ';colorIdx=' num2str(kk)])
   
   rectCoord(1)=min(rectCoordAux(1),rectCoord(1));
   rectCoord(2)=max(rectCoordAux(2),rectCoord(2));
   rectCoord(3)=min(rectCoordAux(3),rectCoord(3));
   rectCoord(4)=max(rectCoordAux(4),rectCoord(4));
end
colormap(jet);
colorbar;
hold(h,'off');

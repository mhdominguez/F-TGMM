%%
%{
%Input:
    coeffs: 1x9 array. Coeffs(7:9)=center of the ellipsoid.
    Coeffs(1:6)=inverse of covariance (symmetric) that defines the ellipsoid.

    kx,ky,kz: scalars to indicate the center of the intersection of XY,YZ,XZ
 
    offsetX,offsetY,offsetZ: scalars. sizeWindowX=1+2*offsetX

    hxy,hxz,hyz: handles to teh axes where you ared isplaying the planes.

    color_: 'g','r','k'

    prymaidLeveleValue: scalar (set to 1 by default)
%}
function handlesDrawEllipse=drawEllipseOrthogonalPlanes(coeffs,kx,ky,kz,offsetX,offsetY,offsetZ,hxy,hxz,hyz,color_,pyramidLevelValue)

scaleSigma=2.0;

%remove this if you want to see contours
%handlesDrawEllipse=drawProbabilityOrthogonalPlanes(coeffs,kx,ky,kz,offsetX,offsetY,offsetZ,hxy,hxz,hyz,color_,pyramidLevelValue);
%return;

%k are the coordinates where plane intersects in x y z
center=coeffs(7:9)+1;%to translate to matlab indexing
handlesDrawEllipse=zeros(1,3);
%intersection of plane and ellipsoid is an ellipse. Equations can be
%written down easily, specially when plane is X=ct or Y=ct or Z=ct;
sigma1=zeros(3);
sigma1(1,1)=coeffs(1);sigma1(2,1)=coeffs(2);sigma1(3,1)=coeffs(3);
sigma1(2,2)=coeffs(4);sigma1(3,2)=coeffs(5);sigma1(3,3)=coeffs(6);
sigma1(1,2)=sigma1(2,1);sigma1(1,3)=sigma1(3,1);sigma1(2,3)=sigma1(3,2);
ll=20;%number of points per ellipse
aa=cos(linspace(0,2.2*pi,ll));
bb=sin(linspace(0,2.2*pi,ll));



switch pyramidLevelValue
    case 1
        
    case 2
        center([1 2])=(center([1 2])+1)/2.0;%+1 is vecause of Matlab indexing in pyramid. 1:2:end->pixel3 goes to pixel2=0.5*(3+1) in next pyramid level
        aux=[2.0 2.0 1.0];
        sigma1=(aux'*aux).*sigma1;
    case 3
        center([1 2])=(center([1 2])+1)/4.0;%in the first two levels only x,y are downsampled
        aux=[4.0 4.0 1.0];
        sigma1=(aux'*aux).*sigma1;
    case 4
        center([1 2])=(center([1 2])+1)/8.0;
        center(3)=(center(3)+1)/2.0;
        aux=[8.0 8.0 2.0];
        sigma1=(aux'*aux).*sigma1;
    case 5
        center([1 2])=(center([1 2])+1)/16.0;
        center(3)=(center(3)+1)/4.0;
        aux=[16.0 16.0 4.0];
        sigma1=(aux'*aux).*sigma1;
    otherwise
end

sigma1=sigma1/scaleSigma;

%YZ intersection
idx=[2 3];
idxN=1;
sigmaX=sigma1(idx,idx);
[V D]=eig(sigmaX);
b=(sigmaX\[sigma1(1,2);sigma1(1,3)])*(kx-center(1));
centerX=center(idx)'-b;

r2=-(kx-center(idxN))^2*sigma1(idxN,idxN)+2*(kx-center(idxN))*(sigma1(idxN,idx(1))*center(idx(1))+sigma1(idxN,idx(2))*center(idx(2)))+1+...
    sigma1(idx(1),idx(1))*b(1)*(b(1)-2*center(idx(1)))+sigma1(idx(2),idx(2))*b(2)*(b(2)-2*center(idx(2)))+...
    2*sigma1(idx(1),idx(2))*(b(1)*b(2)-center(idx(1))*b(2)-center(idx(2))*b(1));
if(r2>0)
    auxPts=V*diag(sqrt(r2./diag(D)))*[aa ;bb];
    auxPts=auxPts'+repmat(centerX',[ll 1]);
    %axes(hyz);
    hold(hyz,'on');
    handlesDrawEllipse(1)=plot(hyz,auxPts(:,2)-offsetZ+1,auxPts(:,1)-offsetY+1,color_,'Linewidth',2);
end
%XY intersection
idx=[1 2];
idxN=3;
sigmaX=sigma1(idx,idx);
[V D]=eig(sigmaX);
b=(sigmaX\[sigma1(1,3);sigma1(2,3)])*(kz-center(3));
centerX=center(idx)'-b;


r2=-(kz-center(idxN))^2*sigma1(idxN,idxN)+2*(kz-center(idxN))*(sigma1(idxN,idx(1))*center(idx(1))+sigma1(idxN,idx(2))*center(idx(2)))+1+...
    sigma1(idx(1),idx(1))*b(1)*(b(1)-2*center(idx(1)))+sigma1(idx(2),idx(2))*b(2)*(b(2)-2*center(idx(2)))+...
    2*sigma1(idx(1),idx(2))*(b(1)*b(2)-center(idx(1))*b(2)-center(idx(2))*b(1));
if(r2>0)
    auxPts=V*diag(sqrt(r2./diag(D)))*[aa ;bb];
    auxPts=auxPts'+repmat(centerX',[ll 1]);
    %axes(hxy);
    hold(hxy,'on');
    handlesDrawEllipse(2)=plot(hxy,auxPts(:,2)-offsetY+1,auxPts(:,1)-offsetX+1,color_,'Linewidth',2);
end

%XZ intersection
idx=[1 3];
idxN=2;
sigmaX=sigma1(idx,idx);
[V D]=eig(sigmaX);
b=(sigmaX\[sigma1(1,2);sigma1(2,3)])*(ky-center(2));
centerX=center(idx)'-b;

r2=-(ky-center(idxN))^2*sigma1(idxN,idxN)+2*(ky-center(idxN))*(sigma1(idxN,idx(1))*center(idx(1))+sigma1(idxN,idx(2))*center(idx(2)))+1+...
    sigma1(idx(1),idx(1))*b(1)*(b(1)-2*center(idx(1)))+sigma1(idx(2),idx(2))*b(2)*(b(2)-2*center(idx(2)))+...
    2*sigma1(idx(1),idx(2))*(b(1)*b(2)-center(idx(1))*b(2)-center(idx(2))*b(1));
if(r2>0)
    auxPts=V*diag(sqrt(r2./diag(D)))*[aa ;bb];
    auxPts=auxPts'+repmat(centerX',[ll 1]);
    %axes(hxz);
    hold(hxz,'on');
    handlesDrawEllipse(3)=plot(hxz,auxPts(:,2)-offsetZ+1,auxPts(:,1)-offsetX+1,color_,'Linewidth',2);
end

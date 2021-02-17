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
function handlesDrawEllipse=drawProbabilityOrthogonalPlanes(coeffs,kx,ky,kz,offsetX,offsetY,offsetZ,hxy,hxz,hyz,color_,pyramidLevelValue)


%color trasnaltions
colorStruct.y=[1 1 0];
colorStruct.m=[1 0 1];
colorStruct.c=[0 1 1];
colorStruct.r=[1 0 0];
colorStruct.g=[0 1 0];
colorStruct.b=[0 0 1];
colorStruct.w=[1 1 1];
colorStruct.k=[0 0 0];


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

wX=kx-offsetX;
wY=ky-offsetY;
wZ=kz-offsetZ;
[XI YI ZI]=ndgrid(kx-wX:kx+wX,ky-wY:ky+wY,kz-wZ:kz+wZ);
[sX sY sZ]=size(XI);
%calculate Gaussian distribution (without normalization)
mu=[XI(:)-center(1) YI(:)-center(2) ZI(:)-center(3)];
boxProb=reshape(exp(-0.5*sum((mu*sigma1).*mu,2)),[sX sY sZ]);
cc=ceil((size(boxProb)+1)/2);
mapC=getfield(colorStruct, color_);

%YZ
hold(hyz, 'on');
colorMask=cat(3, mapC(1)*ones([sY sZ]), mapC(2)*ones([sY sZ]), mapC(3)*ones([sY sZ]));
handlesDrawEllipse(1)=imagesc(colorMask,'Parent',hyz);
set(handlesDrawEllipse(1), 'AlphaData', squeeze(boxProb(cc(1),:,:)));

%XY
hold(hxy, 'on');
colorMask=cat(3, mapC(1)*ones([sX sY]), mapC(2)*ones([sX sY]), mapC(3)*ones([sX sY]));
handlesDrawEllipse(2)=imagesc(colorMask,'Parent',hxy);
set(handlesDrawEllipse(2), 'AlphaData', squeeze(boxProb(:,:,cc(3))));

%XZ
hold(hxz, 'on');
colorMask=cat(3, mapC(1)*ones([sX sZ]), mapC(2)*ones([sX sZ]), mapC(3)*ones([sX sZ]));
handlesDrawEllipse(3)=imagesc(colorMask,'Parent',hxz);
set(handlesDrawEllipse(3), 'AlphaData', squeeze(boxProb(:,cc(2),:)));




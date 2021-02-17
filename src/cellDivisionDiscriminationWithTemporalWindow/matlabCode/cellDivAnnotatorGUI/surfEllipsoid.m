function [X Y Z]=surfEllipsoid(coeffs,numPtsSurf,pyramidLevelValue)

if(isempty(numPtsSurf))
    numPtsSurf=20;
end

if(length(coeffs)~=9)
    error 'Code is only ready for 3D ellipsoids'
end

aux=coeffs;
sigma1(1,1)=aux(1);sigma1(2,1)=aux(2);sigma1(3,1)=aux(3);
sigma1(2,2)=aux(4);sigma1(3,2)=aux(5);sigma1(3,3)=aux(6);
sigma1(1,2)=sigma1(2,1);sigma1(1,3)=sigma1(3,1);sigma1(2,3)=sigma1(3,2);
center=aux(7:9);


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


[V D]=eig(sigma1);
%dead cells sometimes have negative values and crash everyhting
for ii=1:size(D,1)
    if(D(ii,ii)<1e-5) D(ii,ii)=1e-5;
    end
end


[xi yi zi]=ellipsoid(0,0,0,1,1,1,numPtsSurf);
auxPts=V*diag(1.0./sqrt(diag(D)))*[xi(:) yi(:) zi(:)]';
ll=size(auxPts,2);
auxPts=auxPts'+repmat(center,[ll 1]);

X=reshape(auxPts(:,1),[numPtsSurf+1 numPtsSurf+1]);
Y=reshape(auxPts(:,2),[numPtsSurf+1 numPtsSurf+1]);
Z=reshape(auxPts(:,3),[numPtsSurf+1 numPtsSurf+1]);
%============================================================
%draws projection of ellipse into 3 orhtogonal planes
%blob: structure containing all information about blob
%stack: 3D volume with image values
%hXX: handles to 3 axes to plot orthographic planes
function [centerPatch handlesDrawEllipse]=drawOrthogonalPlanes(blob,stack,hxy,hxz,hyz,ww,resetAxes,pyramidLevelValue, offsetCenter)

if(resetAxes>0)
    cla(hxy,'reset');cla(hxz,'reset');cla(hyz,'reset');%it improves speed at least with remote desktop
    
else
    %cla(hxy);cla(hxz);cla(hyz);%it improves speed at least with remote desktop
end


%TODO: padarray to always display the same size window

center=blob(3:5)+1;%to translate to Matlab indexing

switch pyramidLevelValue
    case 1
        
    case 2
        center([1 2])=(center([1 2])+1)/2.0;
        
    case 3
        center([1 2])=(center([1 2])+1)/4.0;%in the first two levels only x,y are downsampled
    case 4
        center([1 2])=(center([1 2])+1)/8.0;
        center(3)=(center(3)+1)/2.0;
    case 5
        center([1 2])=(center([1 2])+1)/16.0;
        center(3)=(center(3)+1)/4.0;
    otherwise
end
center=round(center) + offsetCenter;

imSize=size(stack);
if (isempty(stack) )
    imSize = zeros(1,3);
end

minW=zeros(size(center));
maxW=minW;
for ii=1:3
    minW(ii)=max(center(ii)-ww(ii),1);
    maxW(ii)=min(center(ii)+ww(ii),imSize(ii));
end

prePadding= double(-min([center-ww-1],0));
postPadding=double(max([center+ww]-imSize,0));

patch=stack(minW(1):maxW(1),minW(2):maxW(2),minW(3):maxW(3));

if(sum(prePadding)>0)
    patch=padarray(patch,prePadding,'replicate','pre');
end
if(sum(postPadding)>0)
    patch=padarray(patch,postPadding,'replicate','post');
end

centerPatch=center;

cc=ceil((size(patch)+1)/2);



%axes(hxy);%never use axes: it is definitely much slower
hold(hxy, 'off');
imagesc(squeeze(patch(:,:,cc(3))),'Parent',hxy);colormap gray;
%axes(hxz);
hold(hxz, 'off');
imagesc(squeeze(patch(:,cc(2),:)),'Parent',hxz);colormap gray;
%axes(hyz);
hold(hyz, 'off');
imagesc(squeeze(patch(cc(1),:,:)),'Parent',hyz);colormap gray;

handlesDrawEllipse=drawEllipseOrthogonalPlanes(blob([11:16 3:5]),center(1),center(2),center(3),...
                            center(1)-ww(1),center(2)-ww(2),center(3)-ww(3),hxy,hxz,hyz,'g',pyramidLevelValue);


function [handlesDrawSupervoxels] = drawSupervoxelsBlobSuperimposed(blob,imSize,hxy,hxz,hyz,ww,resetAxes,pyramidLevelValue, offsetCenter, stackIniCache)
%superimposes contours of supervoxels into existing figure
%blob: structure containing all information about blob
%stack: 3D volume with image values
%hXX: handles to 3 axes to plot orthographic planes

global svStructGlobal;

if( sum(imSize) == 0)
    handlesDrawSupervoxels = [];
    return;
end

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
        
        error 'Code is not ready yet'
    case 3
        center([1 2])=(center([1 2])+1)/4.0;%in the first two levels only x,y are downsampled
    case 4
        center([1 2])=(center([1 2])+1)/8.0;
        center(3)=(center(3)+1)/2.0;
        
        error 'Code is not ready yet'
    case 5
        center([1 2])=(center([1 2])+1)/16.0;
        center(3)=(center(3)+1)/4.0;
        
        error 'Code is not ready yet'
    otherwise
        error 'Code is not ready yet'
end
center=round(center) + offsetCenter;

%offset for the drawing box
offsetXYZ = center-ww - 1;
cc = ww + 1;

%draw boundaries of each supervoxel
color_ ={'r','m','c','y'};
colorN = length(color_);

blob_svIdx = blob(19:30);
blob_svIdx = blob_svIdx( blob_svIdx >= 0 );
handlesDrawSupervoxels = -ones(length(blob_svIdx), 3);
for kk = 1:length(blob_svIdx)
    PixelIdxList = svStructGlobal{blob(8) + 1 - stackIniCache + 1}{blob_svIdx(kk)+1};%C-indexing to Matlab indexing
    [sx, sy, sz] = ind2sub(imSize, PixelIdxList+1); %C-indexing to Matlab indexing
    
    
    
    %{
    %drawing a cross on each pixel belonging to the supervoxel (dense display)
    colorAux = [color_{ mod(kk,colorN) + 1 } '+'];
    %XY plane
    pos = find( sz == (cc(3) + offsetXYZ(3)) );
    if( isempty(pos) == false)
        hold(hxy, 'on');
        handlesDrawSupervoxels(kk,1) = plot(hxy, sy(pos) - offsetXYZ(2), sx(pos) - offsetXYZ(1), colorAux);        
    end    
    %XZ plane
    pos = find( sy == (cc(2) + offsetXYZ(2)) );
    if( isempty(pos) == false)
        hold(hxz, 'on');
        handlesDrawSupervoxels(kk,2) = plot(hxz,sz(pos) - offsetXYZ(3), sx(pos) - offsetXYZ(1), colorAux);
    end    
    %YZ plane
    pos = find( sx == (cc(1) + offsetXYZ(1)) );
    if( isempty(pos) == false)
        hold(hyz, 'on');
        handlesDrawSupervoxels(kk,3) = plot(hyz,sz(pos) - offsetXYZ(3), sy(pos) - offsetXYZ(2), colorAux);
    end
    %}
    
    %drawing  just the boundary
    colorAux = [color_{ mod(kk,colorN) + 1} '+'];
    %XY plane
    pos = find( sz == (cc(3) + offsetXYZ(3)) );
    patchBW = false( 2* ww(1) + 1, 2* ww(2) + 1 );
    pos1 = sy(pos) - offsetXYZ(2); pos1 = min(pos1, size(patchBW,1)); pos1 = max(pos1, 1);
    pos2 = sx(pos) - offsetXYZ(1); pos2 = min(pos2, size(patchBW,2)); pos2 = max(pos2, 1);
    
    idx = sub2ind(size(patchBW),pos1, pos2);
    patchBW(idx) = true;
    B = bwboundaries(patchBW,8,'noholes');
    P = [];
    for aa = 1:length(B)
        P = [P;B{aa}];
    end
    if( isempty(pos) == false)
        hold(hxy, 'on');
        handlesDrawSupervoxels(kk,1) = plot(hxy, P(:,1),P(:,2), colorAux);        
    end    
    %XZ plane
    pos = find( sy == (cc(2) + offsetXYZ(2)) );
    patchBW = false( 2* ww(3) + 1, 2* ww(1) + 1 );    
    pos1 = sz(pos) - offsetXYZ(3); pos1 = min(pos1, size(patchBW,1)); pos1 = max(pos1, 1);
    pos2 = sx(pos) - offsetXYZ(1); pos2 = min(pos2, size(patchBW,2)); pos2 = max(pos2, 1);
    
    idx = sub2ind(size(patchBW),pos1, pos2);
    patchBW(idx) = true;
    B = bwboundaries(patchBW,8,'noholes');
    P = [];
    for aa = 1:length(B)
        P = [P;B{aa}];
    end
    if( isempty(pos) == false)
        hold(hxz, 'on');
        handlesDrawSupervoxels(kk,2) = plot(hxz, P(:,1),P(:,2), colorAux);        
    end   
    %YZ plane
    pos = find( sx == (cc(1) + offsetXYZ(1)) );
    patchBW = false( 2* ww(3) + 1, 2* ww(2) + 1 );    
    pos1 = sz(pos) - offsetXYZ(3); pos1 = min(pos1, size(patchBW,1)); pos1 = max(pos1, 1);
    pos2 = sy(pos) - offsetXYZ(2); pos2 = min(pos2, size(patchBW,2)); pos2 = max(pos2, 1);
    
    idx = sub2ind(size(patchBW),pos1, pos2);
    patchBW(idx) = true;
    B = bwboundaries(patchBW,8,'noholes');
    P = [];
    for aa = 1:length(B)
        P = [P;B{aa}];
    end
    if( isempty(pos) == false)
        hold(hyz, 'on');
        handlesDrawSupervoxels(kk,3) = plot(hyz, P(:,1),P(:,2), colorAux);        
    end   
end




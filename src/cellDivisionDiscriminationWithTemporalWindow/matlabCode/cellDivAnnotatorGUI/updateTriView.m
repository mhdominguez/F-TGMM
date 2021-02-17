%======================================================================
%update triview
function handles=updateTriView(handles,frameId,blobId,isDownTheLineage,hObject)

global stackGlobal;
global blobStructGlobal;
global neighGlobal;

ww=get(handles.sliderWindowRadius,'Value');
switch get(handles.popupmenuPyramidLevel,'Value')
    case 1
        ww=round([ww ww ww/handles.anisotropy]);
    case 2
        ww=round([ww ww 2.0*ww/handles.anisotropy]);
    case 3
        ww=round([ww ww ww]);
    case 4
        ww=round([ww ww ww]);
    case 5
        ww=round([ww ww ww]);
    otherwise
end



%update sliders
set(handles.sliderFrame,'Value',frameId);
set(handles.sliderBlob,'Value',blobId);
set(handles.textFrameSlider,'String',num2str(get(handles.sliderFrame,'Value')));
set(handles.textBlobSlider,'String',num2str(get(handles.sliderBlob,'Value')));
set(handles.textObject,'String',['Object: FrameId=' num2str(frameId) ';blobId=' num2str(blobId) ';xyz=[' num2str(blobStructGlobal{frameId+1}(blobId+1,3:5)) ']' ]);

%check if cell is dead
if(blobStructGlobal{frameId+1}(blobId+1,3)<-1e31)
    set(handles.messageText,'String','Current object is dead');
    cla(handles.axes1XY);cla(handles.axes2XZ);cla(handles.axes3YZ);
    cla(handles.axes4XY);cla(handles.axes5XZ);cla(handles.axes6YZ);
    drawnow();
    return;
end

%check if we need to load solution    
pathname = handles.pathLogFile;
pathname = [pathname 'XML_finalResult_lht' filesep];
cacheTGMMxmlSolution(pathname, frameId , handles.frameIni, 0 );%TODO: fix nearest neighbors value from GUI

if(frameId+1<handles.stackIniCache || frameId+1>=handles.stackFinCache)
    %recached set of images
    handles = cacheStackImages(handles,frameId+1,hObject);
    
    %recached supervoxels    
    handles.svFilename = cacheStackSupervoxels(handles, frameId+1, pathname);
        
end

[cc1 handles.drawEllipseAxes123]=drawOrthogonalPlanes(blobStructGlobal{frameId+1}(blobId+1,:),stackGlobal{frameId+1-handles.stackIniCache+1},...
    handles.axes1XY,handles.axes2XZ,handles.axes3YZ,ww,handles.resetAxes,get(handles.popupmenuPyramidLevel,'Value'), handles.offsetCenterPlane123);
handles.resetAxes=max(handles.resetAxes-1,0);
handles.centerTriViewAxes123=cc1;

%superimpose supervoxels
handles.drawSupervoxelAxes123=[];
if( get( handles.checkboxShowBlobSupervoxels, 'Value' ) ~= 0 )
    [handles.drawSupervoxelAxes123] = drawSupervoxelsBlobSuperimposed(blobStructGlobal{frameId+1}(blobId+1,:),size(stackGlobal{frameId+1-handles.stackIniCache+1}),...
                                        handles.axes1XY,handles.axes2XZ,handles.axes3YZ,ww,handles.resetAxes,get(handles.popupmenuPyramidLevel,'Value'), handles.offsetCenterPlane123, handles.stackIniCache);    
end


%check for children
ch=blobStructGlobal{frameId+1}(blobId+1,[17 18]);
if( ch(1) < 0 )
    ch = [];
elseif( ch(2) < 0)%one child
    ch = [frameId + 1, ch(1)];
else%two children
    ch = [frameId + 1, ch(1), frameId+1, ch(2)];
end
if(isempty(ch) )%no children
    set(handles.messageText,'String','Current object does not have chidren');
    cla(handles.axes4XY);cla(handles.axes5XZ);cla(handles.axes6YZ);
    drawnow();
    %check if cell is dead
    
    %=====================extra drawing for dead cells==========================
    %code to display the same region in thenext frame and neighbors, in
    %case we wnat to have an idea of why things dissapear
    bb=blobStructGlobal{frameId+1}(blobId+1,:);
    bb(8)=bb(8)+1;
    cc=drawOrthogonalPlanes(bb,stackGlobal{frameId+2-handles.stackIniCache+1},...
        handles.axes4XY,handles.axes5XZ,handles.axes6YZ,ww,handles.resetAxes,get(handles.popupmenuPyramidLevel,'Value'), handles.offsetCenterPlane456);
    handles.resetAxes=max(handles.resetAxes-1,0);
    handles.centerTriViewAxes456=cc;
    drawEllipseOrthogonalPlanes(blobStructGlobal{frameId+1}(blobId+1,[11:16 3:5]),...
        cc(1),cc(2),cc(3),cc(1)-ww(1),cc(2)-ww(2),cc(3)-ww(3),...
        handles.axes4XY,handles.axes5XZ,handles.axes6YZ,'r',get(handles.popupmenuPyramidLevel,'Value'));%redraw parent
    set(handles.messageText,'String','Current object has no children. Showing next frame for visualization.');
    drawnow();
    
    %draw children of neighbors if check box is on
    
    if(get(handles.checkboxNeighbors,'Value')==true && size( blobStructGlobal,1) >= frameId+2)%make sure I am not in the last frame
        
        xyz = blobStructGlobal{frameId+2}(:,3:5);
        xyz(:,3) = xyz(:,3) * handles.anisotropy;
        xyz( xyz(1)<0,:) = [];
        xyzN = size(xyz,1);
        
        xyzCC=blobStructGlobal{frameId+1}(blobId+1,3:5);
        
        xyzDist=sum((xyz-repmat(xyzCC,[xyzN 1])).^2,2);
        [xyzDist, xyzIdx]=sort(xyzDist,'ascend');
        
        
        for kk=1:min(50,length(xyzIdx)) %numNeighbors we want to display
            
            drawEllipseOrthogonalPlanes(blobStructGlobal{frameId+2}(xyzIdx(kk),[11:16 3:5]),...
                cc1(1),cc1(2),cc1(3),cc1(1)-ww(1),cc1(2)-ww(2),cc1(3)-ww(3),...
                handles.axes4XY,handles.axes5XZ,handles.axes6YZ,'c',get(handles.popupmenuPyramidLevel,'Value'));
            
        end
    end
    %=====================end of extra drawing for dead cells==========================
elseif(blobStructGlobal{ch(1)+1}(ch(2)+1,3)<-1e31)
    set(handles.messageText,'String','Current object children is dead');
    cla(handles.axes4XY);cla(handles.axes5XZ);cla(handles.axes6YZ);
    drawnow();
else
    cc=drawOrthogonalPlanes(blobStructGlobal{ch(1)+1}(ch(2)+1,:),stackGlobal{ch(1)+1-handles.stackIniCache+1},...
        handles.axes4XY,handles.axes5XZ,handles.axes6YZ,ww,handles.resetAxes,get(handles.popupmenuPyramidLevel,'Value'), handles.offsetCenterPlane456);
    handles.resetAxes=max(handles.resetAxes-1,0);
    handles.centerTriViewAxes456=cc;
    drawEllipseOrthogonalPlanes(blobStructGlobal{frameId+1}(blobId+1,[11:16 3:5]),...
        cc(1),cc(2),cc(3),cc(1)-ww(1),cc(2)-ww(2),cc(3)-ww(3),...
        handles.axes4XY,handles.axes5XZ,handles.axes6YZ,'r',get(handles.popupmenuPyramidLevel,'Value'));%redraw parent
    set(handles.messageText,'String','Current object has one chid');
    drawnow();
    if(length(ch)>2)%cell division
        drawEllipseOrthogonalPlanes(blobStructGlobal{ch(3)+1}(ch(4)+1,[11:16 3:5]),...
            cc(1),cc(2),cc(3),cc(1)-ww(1),cc(2)-ww(2),cc(3)-ww(3),...
            handles.axes4XY,handles.axes5XZ,handles.axes6YZ,'g',get(handles.popupmenuPyramidLevel,'Value'));%draw second children
        set(handles.messageText,'String','Current object has two chidren');
        drawnow();
    end
end
%draw neighbors if check box is on
if(get(handles.checkboxNeighbors,'Value')==true)
    
    numNN = round( str2double( get(handles.edit5numNN, 'String')) );
    
    if( numNN > 0 )
        
        nnFrame = neighGlobal{frameId+1};
        if( size(nnFrame,2)  < 2 * numNN )
            %recalculate nearest neighbors within the frame
            xyz = blobStructGlobal{frameId+1}(:,3:5);
            xyz(:,3) = xyz(:,3) * handles.anisotropy;
            idx = knnsearch(xyz,xyz,'k',numNN+1);
            neighGlobal{frameId+1} = frameId * ones( size(idx,1), 2*numNN);
            for jj = 1:size(idx,1)
                neighGlobal{frameId+1}(jj,2:2:end) = idx(jj,2:end)-1;%C-indexing
            end
        end
        
        
        neigh = neighGlobal{frameId+1}(blobId+1,1:2*numNN);
        neigh = reshape(neigh,[2 numNN])';
        pos=find(neigh(:,1)==frameId);
        for kk=1:length(pos)
            drawEllipseOrthogonalPlanes(blobStructGlobal{neigh(pos(kk),1)+1}(neigh(pos(kk),2)+1,[11:16 3:5]),...
                cc1(1),cc1(2),cc1(3),cc1(1)-ww(1),cc1(2)-ww(2),cc1(3)-ww(3),...
                handles.axes1XY,handles.axes2XZ,handles.axes3YZ,'c',get(handles.popupmenuPyramidLevel,'Value'));
        end
        pos=find(neigh(:,1)>frameId);
        for kk=1:length(pos)
            pp=find(neigh(pos(kk),1)==ch(1:2:end) & neigh(pos(kk),2)==ch(2:2:end));
            if(~isempty(pp)) continue;end;
            if(exist('cc','var')==0)
                cc=drawOrthogonalPlanes(blobStructGlobal{neigh(pos(kk),1)+1}(neigh(pos(kk),2)+1,:),stackGlobal{neigh(pos(kk),1)+1-handles.stackIniCache+1},...
                    handles.axes4XY,handles.axes5XZ,handles.axes6YZ,ww,handles.resetAxes,get(handles.popupmenuPyramidLevel,'Value'), handles.offsetCenterPlane456);
                handles.resetAxes=max(handles.resetAxes-1,0);
                drawEllipseOrthogonalPlanes(blobStructGlobal{frameId+1}(blobId+1,[11:16 3:5]),...
                    cc(1),cc(2),cc(3),cc(1)-ww(1),cc(2)-ww(2),cc(3)-ww(3),...
                    handles.axes4XY,handles.axes5XZ,handles.axes6YZ,'r',get(handles.popupmenuPyramidLevel,'Value'));%redraw parent
            end
            drawEllipseOrthogonalPlanes(blobStructGlobal{neigh(pos(kk),1)+1}(neigh(pos(kk),2)+1,[11:16 3:5]),...
                cc(1),cc(2),cc(3),cc(1)-ww(1),cc(2)-ww(2),cc(3)-ww(3),...
                handles.axes4XY,handles.axes5XZ,handles.axes6YZ,'c',get(handles.popupmenuPyramidLevel,'Value'));
        end
        
    end
end




%update dot to indicate lineage tree position
if(~isempty(handles.lineageTreeDot))
    if(isDownTheLineage)
        handles.lineageTreeZlevel=handles.lineageTreeZlevel+1;
    else
        handles.lineageTreeZlevel=handles.lineageTreeZlevel-1;
    end
    handles.lineageTreeDot=showMarkerLineageTree(blobStructGlobal{frameId+1}(blobId+1,:),handles.axesLineageTree,handles.lineageTreeDot,handles.lineageTreeZlevel);
end

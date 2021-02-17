%displays the lineage tree hat contains tree,blobs
%warning: frame and blob use C++ indexing (starts at 0)

%numSolution=0->currentState;
%numSolution>0->solution number
%h,color_ are not mandaotry variables

%output: rectCoord=[xMin xMax yMin yMax];
function rectCoord=displayPartialLineageTree(Zstep,frameOffset,blobStruct,frame,blob,deltaT,numSolution,anisotropy,h,color_,pyramidLevelValue)

%go back deltaT
ff=frame+1-frameOffset;
bb=blob+1;

rectCoord=[1e32 0 1e32 0];%xMin xMax yMin yMax
if(isempty(blobStruct(ff,bb).frame)) return;end




timeT=0;
if(numSolution==0)%current state
    while(blobStruct(ff,bb).currentState.parentIdx(1)<4e9 && time<deltaT)
        ffAux=blobStruct(ff,bb).currentState.parentIdx(1)+1;
        bb=blobStruct(ff,bb).currentState.parentIdx(2)+1;
        ff=ffAux-frameOffset;
        timeT=timeT+1;
    end
else
    while(blobStruct(ff,bb).solutions(numSolution).parentIdx(1)<4e9 && timeT<deltaT)
        ffAux=blobStruct(ff,bb).solutions(numSolution).parentIdx(1)+1;
        bb=blobStruct(ff,bb).solutions(numSolution).parentIdx(2)+1;
        ff=ffAux-frameOffset;
        timeT=timeT+1;
    end
end

%traverse the tree down 
if(exist('h','var')==0)
    h=figure;
%else
%    axes(h);
end

%deltaT-timeT to make sure all lineages are correctly aligned
queue=[ff+frameOffset bb 0+deltaT-timeT  -1 -1 -1];%[frame blob Z-level parent_center]
numFrames=size(blobStruct,1);
while(~isempty(queue))
    ff=queue(1,1)-frameOffset;
    bb=queue(1,2);
    Zlevel=queue(1,3);
    centerPar=queue(1,4:6);
    queue(1,:)=[];
    
    if(Zlevel>2*deltaT) continue;end;%just partial lineage
    
    if(ff>numFrames) break; end;
    
    %add children to queue
    if(numSolution==0)
        ch=blobStruct(ff,bb).currentState.childrenIdx;
    else
        ch=blobStruct(ff,bb).solutions(numSolution).childrenIdx;
    end
    center=blobStruct(ff,bb).surface.coeffs(7:9);
    
    switch pyramidLevelValue
        case 1
            
        case 2
            center([1 2])=(center([1 2])+1)/2.0;%+1 is vecause of Matlab indexing in pyramid. 1:2:end->pixel3 goes to pixel2=0.5*(3+1) in next pyramid level
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
    
    if(center(1)<-1e31)%dead cell: it can not have children
       continue; 
    end
    for ii=1:2:length(ch)
        if(length(ch)>2)%cell division
            queue=[queue;ch(ii)+1 ch(ii+1)+1 Zlevel+1 center];
        else
            queue=[queue;ch(ii)+1 ch(ii+1)+1 Zlevel+1 center];
        end
    end
    
    
    [X Y Z]=surfEllipsoid(blobStruct(ff,bb).surface.coeffs,[],pyramidLevelValue);
    Z=anisotropy*(Z-center(3))+center(3);
    %X=X-center(1);Y=Y-center(2);Z=anisotropy*(Z-center(3));
    %center=[Xstep*Clevel  0 -Zlevel*Zstep];
    %X=X+center(1);Y=Y+center(2);Z=Z+center(3);
    Z=Z-Zlevel*Zstep;
    
    C=ones(size(X));
    if(exist('color_','var'))
        C=color_*C;
    end
    
    %if(ff==frame+1 && bb==blob+1)
    %    C=-C;
    %end
    hold(h,'on');hs=surf(h,X,Y,Z,C);
    set(hs,'FaceAlpha', 0.7);
    
    rectCoord(1)=min(min(X(:)),rectCoord(1));
    rectCoord(2)=max(max(X(:)),rectCoord(2));
    rectCoord(3)=min(min(Y(:)),rectCoord(3));
    rectCoord(4)=max(max(Y(:)),rectCoord(4));
    
    if(sum(centerPar)>0)
        plot3(h,[center(1) centerPar(1)],[center(2) centerPar(2)],[center(3)-Zlevel*Zstep; centerPar(3)-(Zlevel-1)*Zstep;],'Color',[0.5 0.5 0.1],'Linewidth',3);
    end
end
%set the view
view(h,[26 26]);
xlabel(h,'X');
ylabel(h,'Y');
zlabel(h,'Z');
grid(h,'on');
shading flat;

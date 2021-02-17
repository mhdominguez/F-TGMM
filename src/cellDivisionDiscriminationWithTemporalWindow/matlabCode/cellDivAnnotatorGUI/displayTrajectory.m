%displays the lineage tree hat contains tree,blobs
%warning: frame and blob use C++ indexing (starts at 0)

%numSolution=0->currentState;
%numSolution>0->solution number
%h,color_,offset are not mandaotry variables
function displayTrajectory(blobStruct,frame,blob,numSolution,anisotropy,h,color_)



Zstep=30;
Xstep=40;%TODO: plan ahead the width of the tree (i.e. nmber of cell divisions) so we display it appropiately

if(isempty(blobStruct(frame+1,blob+1).frame)) return;end


%go back to the root of the tree
ff=frame+1;
bb=blob+1;

if(numSolution==0)%current state
    while(blobStruct(ff,bb).currentState.parentIdx(1)<4e9)
        ffAux=blobStruct(ff,bb).currentState.parentIdx(1)+1;
        bb=blobStruct(ff,bb).currentState.parentIdx(2)+1;
        ff=ffAux;
    end
else
    while(blobStruct(ff,bb).solutions(numSolution).parentIdx(1)<4e9)
        ffAux=blobStruct(ff,bb).solutions(numSolution).parentIdx(1)+1;
        bb=blobStruct(ff,bb).solutions(numSolution).parentIdx(2)+1;
        ff=ffAux;
    end
end

%traverse the tree down 
if(exist('h','var')==0)
    h=figure;
%else
%    axes(h);
end


queue=[ff bb 0  -1 -1 -1];%[frame blob Z-level parent_center]
numFrames=size(blobStruct,1);
while(~isempty(queue))
    ff=queue(1,1);
    bb=queue(1,2);
    Zlevel=queue(1,3);
    centerPar=queue(1,4:6);
    queue(1,:)=[];
    
    if(ff>numFrames) break; end;
    
    %add children to queue
    if(numSolution==0)
        ch=blobStruct(ff,bb).currentState.childrenIdx;
    else
        ch=blobStruct(ff,bb).solutions(numSolution).childrenIdx;
    end
    center=blobStruct(ff,bb).surface.coeffs(7:9);
    
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
    

    
    %if(ff==frame+1 && bb==blob+1)
    %    C=-C;
    %end
    hold(h,'on');
    
    if(sum(centerPar)>0)
        plot3(h,[center(1) centerPar(1)],[center(2) centerPar(2)],anisotropy*[center(3) centerPar(3)],'-o','Color',color_,'Linewidth',1,'MarkerSize',10,'MarkerEdgeColor','k','MarkerFaceColor',color_);
    end
end
%set the view
view(h,[0 0]);
xlabel(h,'X');
ylabel(h,'Y');
zlabel(h,'Z');
grid(h,'on');
%shading flat;
disp 'REMEMBER TO USE zlim TO ZOOM IN A PARTICULAR AREA OF THE LINEAGE'
%===================================================================================

function handles=saveAnnotationInMemory(handles,hObject,class)

global blobStructGlobal;

blobId=get(handles.sliderBlob,'Value');
frameId=get(handles.sliderFrame,'Value');

if(length(blobStructGlobal) < frameId + 1 || size(blobStructGlobal{frameId+1},1) < blobId+1)
    return;
end;


blob = blobStructGlobal{frameId+1}(blobId+1,:);

%save object as a surface XML object + imFilename + class 
ll=length(handles.objClassifier);
%objAux=blobStructGlobal(frameId+1,blobId+1).surface;
objAux.type = 'Ellipsoid';
objAux.coeffs = [blob(11:16) blob(3:5)];
objAux.id=1;%specific to ellipsoid
objAux.numCoeffs=length(objAux.coeffs);
objAux.covarianceMatrixSize=3;%specific to 3D
objAux.imFilename=[handles.stackFilename{frameId+1}];
objAux.intensity=blob(2);

switch class
    case 0
        objAux.class='nocell';
    case 1
        objAux.class='onecell';
    case 2
        objAux.class='twocell';
    case 12
        objAux.class='halfcell';
    %cases to save mother + daughter1 + daughter2 (for cell division analysis) 
    case 3
        objAux.class='cellDivisionCorrect';%real cell division correctly detected
    case 4
        objAux.class='cellDivisionShouldBe';%cell has not divided when it should
    case 5
        objAux.class='cellDivisionWrong';%cell division that is not a real one (although necessary usually to fix an undersegmentation problem)
    otherwise
        objAux.class='undefined';
end

if(isempty(handles.objClassifier)) 
    handles.objClassifier=objAux;
else
    handles.objClassifier(ll+1)=objAux;
end

%for cases 3,4 and 5 we need to save daughters as well
if( class == 3 || class == 4 || class == 5)
   
    
    ch=blob([17 18]);
    if( ch(1) < 0 )
        ch = [];
    elseif( ch(2) < 0)%one child
        ch = [frameId + 1, ch(1)];
    else%two children
        ch = [frameId + 1, ch(1), frameId+1, ch(2)];
    end
    
    objAux.imFilename=[handles.stackFilename{frameId+2}];
    if(length(ch) > 0)
        %daughter 1
        objAux.type='Ellipsoid';
        objAux.coeffs=blobStructGlobal{ch(1)+1}(ch(2)+1,[11:16 3:5]);
        objAux.intensity=blobStructGlobal{ch(1)+1}(ch(2)+1,2);
    else
        objAux.coeffs = zeros(size( objAux.coeffs ) );
        objAux.intensity = 0;
    end
    ll = ll+1;
    handles.objClassifier(ll+1) = objAux;
    %daughter 2    
    if(length(ch) > 2)
        objAux.type='Ellipsoid';
        objAux.coeffs=blobStructGlobal{ch(3)+1}(ch(4)+1,[11:16 3:5]);
        objAux.intensity=blobStructGlobal{ch(3)+1}(ch(4)+1,2);
    else
        objAux.coeffs = zeros(size( objAux.coeffs ) );
        objAux.intensity = 0;
    end
    ll = ll+1;
    handles.objClassifier(ll+1) = objAux;
end


set(handles.messageText,'String',['Annotation saved as ' objAux.class ' cell.' num2str(ll+1) ' annotations not saved']);

guidata(hObject,handles);

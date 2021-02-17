%call has to be made as follows:
%localLineageDisplay(blobStruct,frameId,blobId,numNN,deltaT,anisotropy,im,frameOffset,pyramidLevelValue)
function varargout = localLineageDisplay(varargin)
% LOCALLINEAGEDISPLAY MATLAB code for localLineageDisplay.fig
%      LOCALLINEAGEDISPLAY, by itself, creates a new LOCALLINEAGEDISPLAY or raises the existing
%      singleton*.
%
%      H = LOCALLINEAGEDISPLAY returns the handle to a new LOCALLINEAGEDISPLAY or the handle to
%      the existing singleton*.
%
%      LOCALLINEAGEDISPLAY('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in LOCALLINEAGEDISPLAY.M with the given input arguments.
%
%      LOCALLINEAGEDISPLAY('Property','Value',...) creates a new LOCALLINEAGEDISPLAY or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before localLineageDisplay_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to localLineageDisplay_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help localLineageDisplay

% Last Modified by GUIDE v2.5 15-Dec-2011 17:22:18

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @localLineageDisplay_OpeningFcn, ...
                   'gui_OutputFcn',  @localLineageDisplay_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before localLineageDisplay is made visible.
function localLineageDisplay_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to localLineageDisplay (see VARARGIN)

% Choose default command line output for localLineageDisplay
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes localLineageDisplay wait for user response (see UIRESUME)
% uiwait(handles.figure1);

%save input variables

%we only need blobStruct at the begining to draw the lineages, so I do not
%need to save it for later callbacks
blobStruct=varargin{1};%It is just a partial version of the whole blobStruct. It should be at least (blobStruct(frameId-deltaT:frameId+deltaT,:)).
handles.frameOffset=varargin{8};% in general frameOffset=frameId-deltaT-1
handles.frameId=varargin{2};
handles.blobId=varargin{3};
handles.numNN=varargin{4};
handles.deltaT=varargin{5};
handles.anisotropy=varargin{6};
handles.im=varargin{7};%image structure containing slices of interest over time points of interest

pyramidLevelValue=varargin{9};


handles.timeT=handles.deltaT+1;%keeps track of the curren time being displayed
handles.sliceZ=round(blobStruct(handles.frameId-handles.frameOffset,handles.blobId).center(3));%keeps track of teh current Z slice being visualized
handles.hSurf=[];

%display lineage trees
[handles.rectCoord handles.Zstep]=neighboringLineageVisualization(handles.frameOffset,blobStruct,handles.frameId,handles.blobId,handles.numNN,handles.deltaT,handles.anisotropy,handles.axesMainFigure,pyramidLevelValue);

%display image
handles=displayImage(handles);

guidata(hObject, handles);


% --- Outputs from this function are returned to the command line.
function varargout = localLineageDisplay_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbuttonUpZ.
function pushbuttonUpZ_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonUpZ (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.sliceZ=handles.sliceZ+1;
handles=displayImage(handles);

guidata(hObject, handles);

% --- Executes on button press in pushbuttonDownZ.
function pushbuttonDownZ_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonDownZ (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.sliceZ=handles.sliceZ-1;
handles=displayImage(handles);

guidata(hObject, handles);

% --- Executes on button press in pushbuttonUpT.
function pushbuttonUpT_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonUpT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.timeT=handles.timeT+1;
handles=displayImage(handles);

guidata(hObject, handles);

% --- Executes on button press in pushbuttonDownT.
function pushbuttonDownT_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonDownT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.timeT=handles.timeT-1;
handles=displayImage(handles);

guidata(hObject, handles);


% --- Executes on button press in checkboxRotate3D.
function checkboxRotate3D_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxRotate3D (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkboxRotate3D

if(get(hObject,'Value')==0)
    rotate3d(handles.axesMainFigure,'off');
else
    rotate3d(handles.axesMainFigure,'on');
end
guidata(hObject, handles);



%------------------------------------------------
%------display image-----------------------
function handles=displayImage(handles)
%check time plane is not out of bounds
numPlanes=length(handles.im);
if(handles.timeT<1)
    handles.timeT=handles.timeT+1;
    return;%we do not need to redraw
elseif(handles.timeT>numPlanes)
    handles.timeT=handles.timeT-1;
    return;
end

%display image
imgzposition=handles.sliceZ-handles.Zstep*(handles.timeT-1);%position in the 3D plot which multiplexes Z-axis between stack slices and different time points
zOffset=handles.im(handles.timeT).rectCoord(5);

%check z-plane is not out of bounds
zIdx=handles.sliceZ-zOffset+1;
zImSize=size(handles.im(handles.timeT).stack,3);
if(zIdx<1)%out of bounds
    handles.sliceZ=handles.sliceZ+1;
    return;%we do not need to redraw
elseif(zIdx>zImSize)%out of bounds
    handles.sliceZ=handles.sliceZ-1;
    return;%we do not need to redraw
end



%draw new image slice
slabSize=str2double(get(handles.editSlabSizeZ,'String'));

if(slabSize<1)
    im=handles.im(handles.timeT).stack(:,:,zIdx);
else
    im=max(handles.im(handles.timeT).stack(:,:,max(zIdx-slabSize,1):min(zIdx+slabSize,zImSize)),[],3);%maximum intensity projection on teh slab
end
handles.hSurf=displayImagePlane(imgzposition,im,handles.axesMainFigure,handles.hSurf,handles.im(handles.timeT).rectCoord);



function editSlabSizeZ_Callback(hObject, eventdata, handles)
% hObject    handle to editSlabSizeZ (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editSlabSizeZ as text
%        str2double(get(hObject,'String')) returns contents of editSlabSizeZ as a double


% --- Executes during object creation, after setting all properties.
function editSlabSizeZ_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editSlabSizeZ (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

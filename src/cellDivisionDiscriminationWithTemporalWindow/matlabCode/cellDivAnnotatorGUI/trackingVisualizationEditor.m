function varargout = trackingVisualizationEditor(varargin)
% TRACKINGVISUALIZATIONEDITOR M-file for trackingVisualizationEditor.fig
%      TRACKINGVISUALIZATIONEDITOR, by itself, creates a new TRACKINGVISUALIZATIONEDITOR or raises the existing
%      singleton*.
%
%      H = TRACKINGVISUALIZATIONEDITOR returns the handle to a new TRACKINGVISUALIZATIONEDITOR or the handle to
%      the existing singleton*.
%
%      TRACKINGVISUALIZATIONEDITOR('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in TRACKINGVISUALIZATIONEDITOR.M with the given input arguments.
%
%      TRACKINGVISUALIZATIONEDITOR('Property','Value',...) creates a new TRACKINGVISUALIZATIONEDITOR or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before trackingVisualizationEditor_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to trackingVisualizationEditor_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help trackingVisualizationEditor

% Last Modified by GUIDE v2.5 10-Mar-2015 09:58:25

global stackGlobal; %we declare it global instead of within the handles since it is read only and it takes a lot of time to copy images around passing them by value
global blobStructGlobal;%cell array. Check readTGMMxmlSolution for more details
global neighGlobal;%cell array. Same structure as blobStrucGlobal in order to store nearest neighbors
global svStructGlobal; %cell array containing the PixelIdxList of all the supervoxels

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @trackingVisualizationEditor_OpeningFcn, ...
                   'gui_OutputFcn',  @trackingVisualizationEditor_OutputFcn, ...
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


% --- Executes just before trackingVisualizationEditor is made visible.
function trackingVisualizationEditor_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to trackingVisualizationEditor (see VARARGIN)

% Choose default command line output for trackingVisualizationEditor
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

%store global variables
handles.maxRAMmemInGB=8;

handles.resetAxes=2;
handles.lineageTreeDot=[];
handles.lineageTreeZlevel=[];
handles.objClassifier=[];
handles.classifierResultsSingleFrame=[];
handles.drawEllipseAxes123=[];
handles.drawEllipseAxes456=[];
handles.stackIniCache=0;
handles.stackFinCache=-1;
handles.stackMaxSizeCache = str2double(get(handles.editMaxCacheSize,'String'));
handles.nextObjectSorted=[-1 -1];%[nextObject currentFrame]. We need current frame to know if we need resort list
handles.nextObjectSortedList=[];%stores the index with order of objects
handles.nextObjectSortedListVal=[];%stores the value with order of objects
handles.centerTriViewAxes123=[];
handles.centerTriViewAxes456=[];

handles.drawSupervoxelAxes123=[];%handles to supervoxels drawing
handles.drawSupervoxelAxes456=[];%handles to supervoxels drawing %TODO: implement supervoxels in bottom triview

handles.offsetCenterPlane123 = [0 0 0];%to be able to move up and down the different planes
handles.offsetCenterPlane456 = [0 0 0];

guidata(hObject, handles);
% UIWAIT makes trackingVisualizationEditor wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = trackingVisualizationEditor_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in load4Dstack.
function load4Dstack_Callback(hObject, eventdata, handles)
% hObject    handle to load4Dstack (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global stackGlobal;
%obtain path to stack
[filename,pathname] = uigetfile('*.txt', 'Select TGMM config file');
if(filename==0)
    return;
end

%extract imgFilePattern
fid = fopen([pathname filename],'r');
while 1
    tline = fgetl(fid);
    if ~ischar(tline), break, end
    
    if( length(tline) > 14 && strcmp(tline(1:14),'imgFilePattern') == 1)
        pos = strfind(tline,'=');
        imgFilePattern = tline(pos(1)+1:end);
    end
    
    if( length(tline) > 11 && strcmp(tline(1:11),'anisotropyZ') == 1)
        pos = strfind(tline,'=');
        handles.anisotropy = str2double( tline(pos(1)+1:end) );
        set(handles.editZanisotropy,'String', num2str(handles.anisotropy));
    end
end
fclose(fid);

%load all the stack

addpath([fileparts( mfilename('fullpath') ) filesep 'readTGMM_XMLoutput' filesep 'generateTiles' ]);

%find first file with extension
for frame = 0:3000
    imgFilename = recoverFilenameFromPattern(imgFilePattern, frame);
    DD = dir([imgFilename '*']);
    if( isempty(DD) == false )
        break;
    end
end

%find number of frames
[~,~,ext] = fileparts(DD(1).name); 

numFrames = 0;
iniFrame = frame;
while( exist([imgFilename ext],'file') > 0 )
    numFrames = numFrames + 1;
    
    frame = frame + 1;
    imgFilename = recoverFilenameFromPattern(imgFilePattern, frame);
end

stackGlobal=cell(handles.stackMaxSizeCache*2+1,1);
handles.stackFilename=cell(numFrames,1);

for ii=iniFrame:frame-1
    imgFilename = recoverFilenameFromPattern(imgFilePattern, ii);
    handles.stackFilename{ii + 1}=[imgFilename ext];
end

rmpath([fileparts( mfilename('fullpath') ) filesep 'readTGMM_XMLoutput' filesep 'generateTiles' ]);

%load images
handles=cacheStackImages(handles,handles.stackMaxSizeCache,hObject);


%load supervoxels
pathname = handles.pathLogFile;
pathname = [pathname 'XML_finalResult_lht' filesep];
handles.svFilename = cacheStackSupervoxels(handles, handles.stackMaxSizeCache, pathname);

%load solution
cacheTGMMxmlSolution(pathname, frame, handles.frameIni, 0 );%TODO: fix nearest neighbors value from GUI

%set string
set(handles.loadedDatasetString,'String',[pathname filename]);



%update the gui
guidata(hObject,handles);
drawnow();

function loadedDatasetString_Callback(hObject, eventdata, handles)
% hObject    handle to loadedDatasetString (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of loadedDatasetString as text
%        str2double(get(hObject,'String')) returns contents of loadedDatasetString as a double

function textFrameSlider_Callback(hObject, eventdata, handles)
% hObject    handle to loadedDatasetString (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of loadedDatasetString as text
%        str2double(get(hObject,'String')) returns contents of loadedDatasetString as a double

global blobStructGlobal;

frameId=round(str2double(get(hObject,'String')));
blobId=get(handles.sliderBlob,'Value');

pathname = handles.pathLogFile;
pathname = [pathname 'XML_finalResult_lht' filesep];
cacheTGMMxmlSolution(pathname, frameId , handles.frameIni, 0 );%TODO: fix nearest neighbors value from GUI

if(frameId>=size(blobStructGlobal,1) || frameId<0 || size(blobStructGlobal{frameId+1},1) < blobId+1)
    set(handles.messageText,'String','Current object does not exist in this solution');
    cla(handles.axes1XY);cla(handles.axes2XZ);cla(handles.axes3YZ);
    cla(handles.axes4XY);cla(handles.axes5XZ);cla(handles.axes6YZ);
    drawnow();
else
    handles=updateTriView(handles,frameId,blobId,true,hObject);

end

set(handles.sliderBlob,'Max',length(blobStructGlobal{frameId+1}-1));
set(handles.sliderFrame,'Value',frameId);
guidata(hObject,handles);


function textBlobSlider_Callback(hObject, eventdata, handles)
% hObject    handle to loadedDatasetString (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of loadedDatasetString as text
%        str2double(get(hObject,'String')) returns contents of loadedDatasetString as a double

global blobStructGlobal;

blobId=round(str2double(get(hObject,'String')));
frameId=get(handles.sliderFrame,'Value');

pathname = handles.pathLogFile;
pathname = [pathname 'XML_finalResult_lht' filesep];
cacheTGMMxmlSolution(pathname, frameId, handles.frameIni,0 );%TODO: fix nearest neighbors value from GUI

if(blobId<0 || blobId>=length(blobStructGlobal{frameId+1}))
    set(handles.messageText,'String','Current object does not exist in this solution');
    cla(handles.axes1XY);cla(handles.axes2XZ);cla(handles.axes3YZ);
    cla(handles.axes4XY);cla(handles.axes5XZ);cla(handles.axes6YZ);
    drawnow();
else
    handles=updateTriView(handles,frameId,blobId,true,hObject);
    

end

%update text display
set(handles.sliderBlob,'Max',length(blobStructGlobal{frameId+1}-1));
set(handles.sliderBlob,'Value',blobId);

guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function loadedDatasetString_CreateFcn(hObject, eventdata, handles)
% hObject    handle to loadedDatasetString (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function messageDisplay_Callback(hObject, eventdata, handles)
% hObject    handle to messageDisplay (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of messageDisplay as text
%        str2double(get(hObject,'String')) returns contents of messageDisplay as a double


% --- Executes during object creation, after setting all properties.
function messageDisplay_CreateFcn(hObject, eventdata, handles)
% hObject    handle to messageDisplay (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function sliderFrame_Callback(hObject, eventdata, handles)
% hObject    handle to sliderFrame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

global blobStructGlobal;

blobId=get(handles.sliderBlob,'Value');
frameId=get(handles.sliderFrame,'Value');

pathname = handles.pathLogFile;
pathname = [pathname 'XML_finalResult_lht' filesep];
cacheTGMMxmlSolution(pathname, frameId , handles.frameIni,0 );%TODO: fix nearest neighbors value from GUI

if(length(blobStructGlobal) < frameId + 1 || size(blobStructGlobal{frameId+1},1) < blobId+1)
    set(handles.messageText,'String','Current object does not exist in this solution');
    cla(handles.axes1XY);cla(handles.axes2XZ);cla(handles.axes3YZ);
    cla(handles.axes4XY);cla(handles.axes5XZ);cla(handles.axes6YZ);
    drawnow();
else
    handles=updateTriView(handles,frameId,blobId,true,hObject);
    
end



set(handles.sliderBlob,'Max',length(blobStructGlobal{frameId+1}-1));
set(handles.textFrameSlider,'String',num2str(get(hObject,'Value')));
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function sliderFrame_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sliderFrame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function sliderBlob_Callback(hObject, eventdata, handles)
% hObject    handle to sliderBlob (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of
%        slider

global blobStructGlobal;

blobId=round(get(handles.sliderBlob,'Value'));
frameId=get(handles.sliderFrame,'Value');

pathname = handles.pathLogFile;
pathname = [pathname 'XML_finalResult_lht' filesep];
cacheTGMMxmlSolution(pathname, frameId , handles.frameIni,0 );%TODO: fix nearest neighbors value from GUI

if(length(blobStructGlobal) < frameId + 1 || size(blobStructGlobal{frameId+1},1) < blobId+1)
    set(handles.messageText,'String','Current object does not exist in this solution');
    cla(handles.axes1XY);cla(handles.axes2XZ);cla(handles.axes3YZ);
    cla(handles.axes4XY);cla(handles.axes5XZ);cla(handles.axes6YZ);
    drawnow();
else
    handles=updateTriView(handles,frameId,blobId,true,hObject);      
end

%update text display
if(isfield(get(hObject),'Value')) %sometime it is called from keyPressFcn and it does not have a value
    set(handles.sliderBlob,'Max',length(blobStructGlobal{frameId+1}-1));
    set(handles.textBlobSlider,'String',num2str(get(hObject,'Value')));
end

%sprintf('Score = %g\n',blobStructGlobal{frameId+1}(blobId+1,9))

guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function sliderBlob_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sliderBlob (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function sliderLabel_Callback(hObject, eventdata, handles)
% hObject    handle to sliderLabel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

set(handles.textLabelSlider,'String',num2str(get(hObject,'Value')));
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function sliderLabel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sliderLabel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in pushButtonNextChildLeft.
function pushButtonNextChildLeft_Callback(hObject, eventdata, handles)
% hObject    handle to pushButtonNextChildLeft (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global blobStructGlobal;

blobId=get(handles.sliderBlob,'Value');
frameId=get(handles.sliderFrame,'Value');

pathname = handles.pathLogFile;
pathname = [pathname 'XML_finalResult_lht' filesep];
cacheTGMMxmlSolution(pathname, frameId , handles.frameIni,0 );%TODO: fix nearest neighbors value from GUI

ch=blobStructGlobal{frameId+1}(blobId+1,[17 18]);
if( ch(1) < 0 )
    ch = [];
elseif( ch(2) < 0)%one child
    ch = [frameId + 1, ch(1)];
else%two children
    ch = [frameId + 1, ch(1), frameId+1, ch(2)];
end

if(length(ch)==2)
    frameId=ch(1);
    blobId=ch(2);  
elseif(length(ch)==4)
    frameId=ch(1);
    blobId=ch(2);  
end

pathname = handles.pathLogFile;
pathname = [pathname 'XML_finalResult_lht' filesep];
cacheTGMMxmlSolution(pathname, frameId , handles.frameIni,0 );%TODO: fix nearest neighbors value from GUI


if(isempty(ch))
    set(handles.messageText,'String','Current object has no children');
else
   handles=updateTriView(handles,frameId,blobId,true,hObject);
end



guidata(hObject,handles);


% --- Executes on button press in pushbuttonNextChildRight.
function pushbuttonNextChildRight_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonNextChildRight (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global blobStructGlobal;

blobId=get(handles.sliderBlob,'Value');
frameId=get(handles.sliderFrame,'Value');

pathname = handles.pathLogFile;
pathname = [pathname 'XML_finalResult_lht' filesep];
cacheTGMMxmlSolution(pathname, frameId, handles.frameIni,0 );%TODO: fix nearest neighbors value from GUI

ch=blobStructGlobal{frameId+1}(blobId+1,[17 18]);
if( ch(1) < 0 )
    ch = [];
elseif( ch(2) < 0)%one child
    ch = [frameId + 1, ch(1)];
else%two children
    ch = [frameId + 1, ch(1), frameId+1, ch(2)];
end

if(length(ch)==2)
    frameId=ch(1);
    blobId=ch(2);  
elseif(length(ch)==4)
    frameId=ch(3);
    blobId=ch(4);  
end

pathname = handles.pathLogFile;
pathname = [pathname 'XML_finalResult_lht' filesep];
cacheTGMMxmlSolution(pathname, frameId , handles.frameIni,0 );%TODO: fix nearest neighbors value from GUI


if(isempty(ch))
    set(handles.messageText,'String','Current object has no children');
else
    handles=updateTriView(handles,frameId,blobId,true,hObject);
end



guidata(hObject,handles);


% --- Executes on button press in pushbuttonParent.
function pushbuttonParent_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonParent (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global blobStructGlobal;

blobId=get(handles.sliderBlob,'Value');
frameId=get(handles.sliderFrame,'Value');

par = [frameId-1, blobStructGlobal{frameId + 1}(blobId+1,7)];
if( par(1) < 0 || par(2) < 0 )
    par(1) = 4294967295;
end


frameId=par(1);
blobId=par(2);    

if(par(1)>4e9)
    set(handles.messageText,'String','Current object has no parent');
else
    pathname = handles.pathLogFile;
    pathname = [pathname 'XML_finalResult_lht' filesep];
    cacheTGMMxmlSolution(pathname, frameId , handles.frameIni,0 );%TODO: fix nearest neighbors value from GUI

    handles=updateTriView(handles,frameId,blobId,false,hObject);
end

guidata(hObject,handles);

% --- Executes on button press in pushbuttonLoadSolution.
function pushbuttonLoadSolution_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonLoadSolution (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


global stackGlobal;
%obtain path to stack
[filename,pathname] = uigetfile('experimentLog*.txt', 'Select TGMM log output file');
if(filename==0)
    return;
end

handles.pathLogFile = pathname;%so we can use it as a reference later to save annotations

%extract imgFilePattern and anisotropyZ
fid = fopen([pathname filename],'r');
while 1
    tline = fgetl(fid);
    if ~ischar(tline), break, end
    
    if( length(tline) > 14 && strcmp(tline(1:14),'imgFilePattern') == 1)
        pos = strfind(tline,'=');
        handles.imgFilePattern = tline(pos(1)+1:end);
    end
    
    if( length(tline) > 11 && strcmp(tline(1:11),'anisotropyZ') == 1)
        pos = strfind(tline,'=');
        handles.anisotropy = str2double( tline(pos(1)+1:end) );
        set(handles.editZanisotropy,'String', num2str(handles.anisotropy));
    end    
    
end
fclose(fid);


%calculate number of frames
pathnameXML = [pathname 'XML_finalResult_lht' filesep];
D = dir([pathnameXML 'GMEMfinalResult_frame*.xml']);
frameIni = str2double(D(1).name(end-7:end-4));
frameEnd = str2double(D(end).name(end-7:end-4));
numFrames = frameEnd - frameIni + 1;

handles.frameIni = frameIni;
guidata(hObject,handles);

%load all the image stacks into cache
%find first file with extension
for frame = frameIni:frameEnd
    imgFilename = recoverFilenameFromPattern(handles.imgFilePattern, frame);
    DD = dir([imgFilename '*']);
    if( isempty(DD) == false )
        break;
    end
end

%find number of frames
[~,~,ext] = fileparts(DD(1).name); 


stackGlobal=cell(handles.stackMaxSizeCache*2+1,1);
handles.stackFilename=cell(numFrames + frameIni + 1,1);

for ii=frameIni:frameEnd+1
    imgFilename = recoverFilenameFromPattern(handles.imgFilePattern, ii);
    handles.stackFilename{ii + 1}=[imgFilename ext];
end

tic;
%load images
handles=cacheStackImages(handles,handles.stackMaxSizeCache + handles.frameIni + 1,hObject);
tt = toc;
disp(['Took ' num2str(tt) ' secs to load ' num2str(sum(~cellfun(@isempty,stackGlobal))) ' images']);
%set string
set(handles.loadedDatasetString,'String',[pathname filename]);

%update the gui
guidata(hObject,handles);
drawnow();

%load specific solution

global blobStructGlobal;
global neighGlobal;
global svStructGlobal;
%global stackGlobal;

% % %[filename,pathname] = uigetfile('*.xml', 'Select solution');
% % 
% % % % if(filename==0)
% % % %     return;
% % % % end
% % guidata(hObject,handles);
% % drawnow();


set(handles.messageText,'String','Loading XML solutiona from TGMM...');
%update the gui
guidata(hObject,handles);
drawnow();


neighGlobal = cell(frameEnd + 1,1);

tic;
blobStructGlobal = cell( frameEnd + 1,1);
cacheTGMMxmlSolution(pathnameXML,frame, frameIni, 0 );%TODO: fix nearest neighbors value from GUI

tt = toc;
disp(['Took ' num2str(tt) ' secs to load 1 xml file']);


%-------------------------------------------
set(handles.messageText,'String','Loading supervoxels...');
%update the gui
guidata(hObject,handles);
drawnow();
tic;
svStructGlobal=cell(handles.stackMaxSizeCache*2+1,1);
handles.svFilename = cacheStackSupervoxels(handles, handles.stackMaxSizeCache  + handles.frameIni + 1, pathnameXML);

tt = toc;
disp(['Took ' num2str(tt) ' secs to load supervoxels from ' num2str(sum(~cellfun(@isempty,svStructGlobal))) ' images']);
%-------------------------------------------------------
%update scrolling bars
%numFrames=size(blobStructGlobal,1);
ll = cellfun(@length,blobStructGlobal);
maxBeads = max(ll);



set(handles.sliderFrame,'Max',frameEnd);
set(handles.textFrameSlider,'String',num2str(frameIni));
set(handles.sliderFrame,'Value',frameIni);
set(handles.sliderFrame,'Min',frameIni);
set(handles.sliderBlob,'Min',0);
set(handles.sliderBlob,'Max',maxBeads-1);
set(handles.textBlobSlider,'String','0');
set(handles.sliderBlob,'Value',0);
set(handles.sliderBlob,'SliderStep',[1.0/(maxBeads-1) 10.0/(maxBeads-1)]);

if(numFrames==1)
    set(handles.sliderFrame,'SliderStep',[0.1 0.1]);
else
    set(handles.sliderFrame,'SliderStep',[1.0/(numFrames-1) 10.0/(numFrames-1)]);
end





set(handles.messageText,'String','Solution loaded succesfully');

guidata(hObject,handles);

%=============================================================================================================
%%
% --- Executes on slider movement.
function sliderWindowRadius_Callback(hObject, eventdata, handles)
% hObject    handle to sliderWindowRadius (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
handles.resetAxes=2;
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function sliderWindowRadius_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sliderWindowRadius (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function sliderAzimuthal_Callback(hObject, eventdata, handles)
% hObject    handle to sliderAzimuthal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

view(handles.axesLineageTree,[get(handles.sliderAzimuthal,'Value') get(handles.sliderElevation,'Value')]);
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function sliderAzimuthal_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sliderAzimuthal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function sliderElevation_Callback(hObject, eventdata, handles)
% hObject    handle to sliderElevation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
view(handles.axesLineageTree,[get(handles.sliderAzimuthal,'Value') get(handles.sliderElevation,'Value')]);
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function sliderElevation_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sliderElevation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in checkboxNeighbors.
function checkboxNeighbors_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxNeighbors (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkboxNeighbors
if(get(hObject,'Value')==true)%option checked
    sliderFrame_Callback(hObject, eventdata, handles);%redo scene with neighbors
else
    %redo scene without neighbors
    sliderFrame_Callback(hObject, eventdata, handles);
end


% --- Executes on slider movement.
function sliderZmax_Callback(hObject, eventdata, handles)
% hObject    handle to sliderZmax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
zlim(handles.axesLineageTree,[get(handles.sliderZmin,'Value') get(handles.sliderZmax,'Value') ]);
zl=zlim(handles.axesLineageTree);
set(handles.sliderZcenter,'Value',sum(zl)/2);
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function sliderZmax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sliderZmax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function sliderZmin_Callback(hObject, eventdata, handles)
% hObject    handle to sliderZmin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of
%        slider
zlim(handles.axesLineageTree,[get(handles.sliderZmin,'Value') get(handles.sliderZmax,'Value') ]);
zl=zlim(handles.axesLineageTree);
set(handles.sliderZcenter,'Value',sum(zl)/2);
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function sliderZmin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sliderZmin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function sliderZcenter_Callback(hObject, eventdata, handles)
% hObject    handle to sliderZcenter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
zl=zlim(handles.axesLineageTree);
zl=(zl(2)-zl(1))/2;
zlim(handles.axesLineageTree,[get(hObject,'Value')-zl get(hObject,'Value')+zl ]);

%update elements
zl=zlim(handles.axesLineageTree);
set(handles.sliderZmin,'Value',zl(1));
set(handles.sliderZmax,'Value',zl(2));
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function sliderZcenter_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sliderZcenter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


%%


% --- Executes on button press in pushbuttonNextObject.
function pushbuttonNextObject_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonNextObject (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global blobStructGlobal;


frameId=get(handles.sliderFrame,'Value');

maxBeads=length(blobStructGlobal{frameId+1}); 
%set(handles.sliderBlob,'SliderStep',[1.0/(maxBeads-1) 100.0/(maxBeads-1)]);%activate to change slider

%blobId=get(handles.sliderBlob,'Value')+1;
blobId=randi(maxBeads)-1;%at random to validate results

if(~isempty(handles.classifierResultsSingleFrame))
    
    while(handles.classifierResultsSingleFrame(blobId+1)<-30)%change this to adjust Fx threshold
        blobId=randi(maxBeads)-1;
    end
    display(['Using classifier results to select blob. Fx=' num2str(handles.classifierResultsSingleFrame(blobId+1))]);
else
    while( size(blobStructGlobal{frameId+1},1) < blobId+1 || blobStructGlobal{frameId+1}(blobId+1,3)<-1e31)%to avoid displaying dead cells
        blobId=randi(maxBeads)-1;
    end
end

set(handles.sliderBlob,'Value',blobId);
guidata(hObject,handles);
sliderBlob_Callback(hObject, eventdata, handles);


set(handles.textBlobSlider,'String',num2str(blobId));
guidata(hObject,handles);


% --- Executes on button press in pushbuttonClassifier0Cell.
function pushbuttonClassifier0Cell_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonClassifier0Cell (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles=saveAnnotationInMemory(handles,hObject,0);
guidata(hObject,handles);


% --- Executes on button press in pushbuttonClassifier1Cell.
function pushbuttonClassifier1Cell_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonClassifier1Cell (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles=saveAnnotationInMemory(handles,hObject,1);
guidata(hObject,handles);

% --- Executes on button press in pushbuttonClassifier2Cell.
function pushbuttonClassifier2Cell_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonClassifier2Cell (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles=saveAnnotationInMemory(handles,hObject,2);
guidata(hObject,handles);

% --- Executes on button press in pushbuttonClassifier1_2Cell.
function pushbuttonClassifier1_2Cell_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonClassifier1_2Cell (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles=saveAnnotationInMemory(handles,hObject,12);

guidata(hObject,handles);
% --- Executes on button press in pushbuttonClassifierSave2Xml.
function pushbuttonClassifierSave2Xml_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonClassifierSave2Xml (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


if(~isempty(handles.objClassifier))
    %mkdir if necessary
    pathAnn = [handles.pathLogFile 'annForCellDivDiscrWithTempWin' filesep];
    
    if( exist(pathAnn,'dir') == 0 )
        mkdir(pathAnn);
    end
    
    %save all the annotations as an XML surface object
    %no header or footer to be able to combine multiple savings with CAT
    %easily        
    fileOut = [pathAnn 'classifierAnnotations_' datestr(now,30) '.xml'];
    fid=fopen(fileOut,'w');
    %write header
    fprintf(fid,'<?xml version="1.0" encoding="UTF-8"?>\n<document>\n');
    %write body
    if( isfield(handles.objClassifier, 'svIdx' ) == 0 )
        writeSurfaceObjectToXML(handles.objClassifier,fid);%we do not have supervoxel information
    else
        writeSurfaceObjectToXMLsupervoxel(handles.objClassifier,fid);
    end
    %write footer
    fprintf(fid,'</document>\n');
    fclose(fid);
    %output message in message center text box
    set(handles.messageText,'String',[ num2str(length(handles.objClassifier)) ' annotations written successfully to file ' fileOut]);
else
    set(handles.messageText,'String','No annotations to be saved at this point');
end

handles.objClassifier=[];
guidata(hObject,handles);


% --- Executes on key press with focus on figure1 and none of its controls.
function figure1_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  structure with the following fields (see FIGURE)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)

switch eventdata.Key
    
    case 'n'
        pushbuttonNextObject_Callback(hObject, eventdata, handles); %select next object in this frame at random
    case 'f'
        %pushbuttonNextSorted_Callback(hObject, eventdata, handles); %select next sorted object in this frame at random
        handles = pushbuttonSorted(hObject,eventdata,handles,true);
        guidata(hObject,handles);
        
    case 'd'
        %pushbuttonPrevSorted_Callback(hObject, eventdata, handles); %select next sorted object in this frame at random
        handles = pushbuttonSorted(hObject,eventdata,handles,false);
        guidata(hObject,handles);
    case {'0','numpad0'}
        handles=saveAnnotationInMemory(handles,hObject,0);
        guidata(hObject,handles);
    case {'1','numpad1'}
        handles=saveAnnotationInMemory(handles,hObject,1);
        guidata(hObject,handles);
    case {'2','numpad2'}
        handles=saveAnnotationInMemory(handles,hObject,2);
        guidata(hObject,handles);
    case 'h' %classify cell as oversegmented (half a cell)
        handles=saveAnnotationInMemory(handles,hObject,12);
        guidata(hObject,handles);        
    case 'l'
        pushButtonNextChildLeft_Callback(hObject, eventdata, handles);%left child
    case 'r'
        pushbuttonNextChildRight_Callback(hObject, eventdata, handles);%right child
    case 'p'
        pushbuttonParent_Callback(hObject, eventdata, handles);%right child        
    case 'w'
        pushbuttonWrongCellDivision_Callback(hObject, eventdata, handles);
    case 'c'
        pushbuttonCorrectCellDivision_Callback(hObject, eventdata, handles);
    case 's'
        pushbuttonShouldBeCellDivision_Callback(hObject, eventdata, handles);
    case 'rightarrow'
        blobId = get(handles.sliderBlob,'Value');
        set(handles.sliderBlob,'Value',blobId + 1);
        guidata(hObject,handles);
        sliderBlob_Callback(hObject, eventdata, handles);
    case 'leftarrow'
        blobId = get(handles.sliderBlob,'Value');
        set(handles.sliderBlob,'Value',blobId - 1);
        guidata(hObject,handles);
        sliderBlob_Callback(hObject, eventdata, handles);
    case 'pagedown' % move to the next blob multiple of 10
        blobId = get(handles.sliderBlob,'Value');
        set(handles.sliderBlob,'Value', ceil (blobId / 10.0) * 10 );
        guidata(hObject,handles);
        sliderBlob_Callback(hObject, eventdata, handles);
    case 'pageup' % move to the next blob multiple of 10
        blobId = get(handles.sliderBlob,'Value');
        set(handles.sliderBlob,'Value', floor(blobId / 10.0) * 10 );
        guidata(hObject,handles);
        sliderBlob_Callback(hObject, eventdata, handles);
end

%handles.prevKey=eventdata.Key;
%guidata(hObject,handles); //we don't always have to do it (only when we do
%not call callbacks)


% --- Executes on button press in pushbuttonClassifierUndo.
function pushbuttonClassifierUndo_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonClassifierUndo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


if(~isempty(handles.objClassifier))
    handles.objClassifier(length(handles.objClassifier))=[];
    set(handles.messageText,'String',['Last annotation removed.' num2str(length(handles.objClassifier)) ' annotations not saved']);
else
    set(handles.messageText,'String',['No annotations to undo']);
end


guidata(hObject,handles);


% --- Executes on button press in pushbuttonLoadClassifierResults.
function pushbuttonLoadClassifierResults_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonLoadClassifierResults (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename,pathname] = uigetfile('*.txt', 'Select one element of classifier results');
if(filename==0)
    return;
end
%the results stack
handles.classifierResultsSingleFrame=load([pathname filename]);

set(handles.messageText,'String','Loaded classifier results successfully');
%update the gui
guidata(hObject,handles);


% --- Executes on button press in pushbuttonEditMerge.
function pushbuttonEditMerge_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonEditMerge (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%update the gui
guidata(hObject,handles);


% --- Executes on button press in pushbuttonEditSplit.
function pushbuttonEditSplit_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonEditSplit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


%update the gui
guidata(hObject,handles);

% --- Executes on button press in pushbuttonEditDelete.
function pushbuttonEditDelete_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonEditDelete (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global blobStructGlobal;

blobId=get(handles.sliderBlob,'Value');
frameId=get(handles.sliderFrame,'Value');

if(isempty(blobStructGlobal(frameId+1,blobId+1).frame) || blobStructGlobal(frameId+1,blobId+1).center(1)<-1e31)
    return;
end

%kill cell
blobStructGlobal(frameId+1,blobId+1).center=repmat(-1e32,size(blobStructGlobal(frameId+1,blobId+1).center));

par=blobStructGlobal(frameId+1,blobId+1).solutions(1).parentIdx;
ch=blobStructGlobal(frameId+1,blobId+1).solutions(1).childrenIdx;

blobStructGlobal(frameId+1,blobId+1).solutions(1).parentIdx=[uint32(2^32-1) 0];
blobStructGlobal(frameId+1,blobId+1).solutions(1).childrenIdx=[];

%break tree ties
for ii=1:2:length(ch)
    blobStructGlobal(ch(ii)+1,ch(ii+1)+1).solutions(1).parentIdx=[uint32(2^32-1) 0];
end

if(par(1)<4e9)
    ch=blobStructGlobal(par(1)+1,par(2)+1).solutions(1).childrenIdx;
    for ii=1:2:length(ch)
        if(ch(ii)==frameId && ch(ii+1)==blobId)
            blobStructGlobal(par(1)+1,par(2)+1).solutions(1).childrenIdx([ii ii+1])=[];
            break;
        end
    end
end

set(handles.messageText,'String','Cell deleted');
%update the gui
guidata(hObject,handles);

% --- Executes on button press in pushbuttonEditAdd.
function pushbuttonEditAdd_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonEditAdd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


%update the gui
guidata(hObject,handles);

% --- Executes on slider movement.
function sliderEditScale_Callback(hObject, eventdata, handles)
% hObject    handle to sliderEditScale (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

global blobStructGlobal;

val=get(hObject,'Value');

if(handles.sliderEditScaleVal<val) %scale up
    sc=1.3;
else%scale down
    sc=1.0/1.3;
end
handles.sliderEditScaleVal=val;%save value for later


%scale blob
blobId=get(handles.sliderBlob,'Value');
frameId=get(handles.sliderFrame,'Value');

if(isempty(blobStructGlobal(frameId+1,blobId+1).frame))
    set(handles.messageText,'String','Current object does not exist in this solution');
    return;
end
blobStructGlobal(frameId+1,blobId+1).surface.coeffs(1:6)=blobStructGlobal(frameId+1,blobId+1).surface.coeffs(1:6)/sc;


%delete current display of blob: we don't need to reload image
if(~isempty(handles.drawEllipseAxes123))
    
    for ii=1:length(handles.drawEllipseAxes123)
        if(handles.drawEllipseAxes123(ii)>0)%if object dissapears it returns 0
        delete(handles.drawEllipseAxes123(ii))
        end
    end
    handles.drawEllipseAxes123=[];
end

%draw new blob
ww=get(handles.sliderWindowRadius,'Value');
switch get(handles.popupmenuPyramidLevel,'Value')
    case 1
        ww=round([ww ww ww/handles.anisotropy]);
    case 2
        
        ww=round([ww ww 2.0*ww/handles.anisotropy]);
    case 3
        
    case 4
        
    case 5
        
    otherwise
end
cc=handles.centerTriViewAxes123;
handles.drawEllipseAxes123=drawEllipseOrthogonalPlanes(blobStructGlobal(frameId+1,blobId+1).surface.coeffs,...
    cc(1),cc(2),cc(3),cc(1)-ww(1),cc(2)-ww(2),cc(3)-ww(3),...
    handles.axes1XY,handles.axes2XZ,handles.axes3YZ,'g',get(handles.popupmenuPyramidLevel,'Value'));%redraw parent

drawnow();

%update the gui
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function sliderEditScale_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sliderEditScale (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

handles.sliderEditScaleVal=get(hObject,'Value');%save value for later


%update the gui
guidata(hObject,handles);

% --- Executes on slider movement.
function sliderEditScaleX_Callback(hObject, eventdata, handles)
% hObject    handle to sliderEditScaleX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

global blobStructGlobal;

val=get(hObject,'Value');

if(handles.sliderEditScaleXVal<val) %scale up
    sc=1.3;
else%scale down
    sc=1.0/1.3;
end
handles.sliderEditScaleXVal=val;%save value for later


%scale blob
blobId=get(handles.sliderBlob,'Value');
frameId=get(handles.sliderFrame,'Value');

if(isempty(blobStructGlobal(frameId+1,blobId+1).frame))
    set(handles.messageText,'String','Current object does not exist in this solution');
    return;
end

coeffs=blobStructGlobal(frameId+1,blobId+1).surface.coeffs;
sigma1=zeros(3);
sigma1(1,1)=coeffs(1);sigma1(2,1)=coeffs(2);sigma1(3,1)=coeffs(3);
sigma1(2,2)=coeffs(4);sigma1(3,2)=coeffs(5);sigma1(3,3)=coeffs(6);
sigma1(1,2)=sigma1(2,1);sigma1(1,3)=sigma1(3,1);sigma1(2,3)=sigma1(3,2);

[U S V]=svd(sigma1,'econ');
S(1,1)=S(1,1)/sc;
sigma1=U*S*V';

blobStructGlobal(frameId+1,blobId+1).surface.coeffs(1:6)=[sigma1(1,1) sigma1(1,2) sigma1(1,3) sigma1(2,2) sigma1(2,3) sigma1(3,3)];


%delete current display of blob: we don't need to reload image
if(~isempty(handles.drawEllipseAxes123))
    
    for ii=1:length(handles.drawEllipseAxes123)
        if(handles.drawEllipseAxes123(ii)>0)%if object dissapears it returns 0
        delete(handles.drawEllipseAxes123(ii))
        end
    end
    handles.drawEllipseAxes123=[];
end

%draw new blob
ww=get(handles.sliderWindowRadius,'Value');
switch get(handles.popupmenuPyramidLevel,'Value')
    case 1
        ww=round([ww ww ww/handles.anisotropy]);
    case 2
        
        ww=round([ww ww 2.0*ww/handles.anisotropy]);
    case 3
        
    case 4
        
    case 5
        
    otherwise
end
cc=handles.centerTriViewAxes123;
handles.drawEllipseAxes123=drawEllipseOrthogonalPlanes(blobStructGlobal(frameId+1,blobId+1).surface.coeffs,...
    cc(1),cc(2),cc(3),cc(1)-ww(1),cc(2)-ww(2),cc(3)-ww(3),...
    handles.axes1XY,handles.axes2XZ,handles.axes3YZ,'g',get(handles.popupmenuPyramidLevel,'Value'));%redraw parent
drawnow();

%update the gui
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function sliderEditScaleX_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sliderEditScaleX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

handles.sliderEditScaleXVal=get(hObject,'Value');%save value for later
%update the gui
guidata(hObject,handles);

% --- Executes on slider movement.
function sliderEditScaleY_Callback(hObject, eventdata, handles)
% hObject    handle to sliderEditScaleY (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


global blobStructGlobal;

val=get(hObject,'Value');

if(handles.sliderEditScaleYVal<val) %scale up
    sc=1.3;
else%scale down
    sc=1.0/1.3;
end
handles.sliderEditScaleYVal=val;%save value for later


%scale blob
blobId=get(handles.sliderBlob,'Value');
frameId=get(handles.sliderFrame,'Value');

if(isempty(blobStructGlobal(frameId+1,blobId+1).frame))
    set(handles.messageText,'String','Current object does not exist in this solution');
    return;
end

coeffs=blobStructGlobal(frameId+1,blobId+1).surface.coeffs;
sigma1=zeros(3);
sigma1(1,1)=coeffs(1);sigma1(2,1)=coeffs(2);sigma1(3,1)=coeffs(3);
sigma1(2,2)=coeffs(4);sigma1(3,2)=coeffs(5);sigma1(3,3)=coeffs(6);
sigma1(1,2)=sigma1(2,1);sigma1(1,3)=sigma1(3,1);sigma1(2,3)=sigma1(3,2);

[U S V]=svd(sigma1,'econ');
S(2,2)=S(2,2)/sc;
sigma1=U*S*V';

blobStructGlobal(frameId+1,blobId+1).surface.coeffs(1:6)=[sigma1(1,1) sigma1(1,2) sigma1(1,3) sigma1(2,2) sigma1(2,3) sigma1(3,3)];


%delete current display of blob: we don't need to reload image
if(~isempty(handles.drawEllipseAxes123))
    
    for ii=1:length(handles.drawEllipseAxes123)
        if(handles.drawEllipseAxes123(ii)>0)%if object dissapears it returns 0           
            delete(handles.drawEllipseAxes123(ii))
        end
    end
    handles.drawEllipseAxes123=[];
end

%draw new blob
ww=get(handles.sliderWindowRadius,'Value');
switch get(handles.popupmenuPyramidLevel,'Value')
    case 1
        ww=round([ww ww ww/handles.anisotropy]);
    case 2
        
        ww=round([ww ww 2.0*ww/handles.anisotropy]);
    case 3
        
    case 4
        
    case 5
        
    otherwise
end
cc=handles.centerTriViewAxes123;
handles.drawEllipseAxes123=drawEllipseOrthogonalPlanes(blobStructGlobal(frameId+1,blobId+1).surface.coeffs,...
    cc(1),cc(2),cc(3),cc(1)-ww(1),cc(2)-ww(2),cc(3)-ww(3),...
    handles.axes1XY,handles.axes2XZ,handles.axes3YZ,'g',get(handles.popupmenuPyramidLevel,'Value'));%redraw parent
drawnow();

%update the gui
guidata(hObject,handles);

%update the gui
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function sliderEditScaleY_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sliderEditScaleY (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

handles.sliderEditScaleYVal=get(hObject,'Value');%save value for later
%update the gui
guidata(hObject,handles);

% --- Executes on slider movement.
function sliderEditScaleZ_Callback(hObject, eventdata, handles)
% hObject    handle to sliderEditScaleZ (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

global blobStructGlobal;

val=get(hObject,'Value');

if(handles.sliderEditScaleZVal<val) %scale up
    sc=1.3;
else%scale down
    sc=1.0/1.3;
end
handles.sliderEditScaleZVal=val;%save value for later


%scale blob
blobId=get(handles.sliderBlob,'Value');
frameId=get(handles.sliderFrame,'Value');

if(isempty(blobStructGlobal(frameId+1,blobId+1).frame))
    set(handles.messageText,'String','Current object does not exist in this solution');
    return;
end

coeffs=blobStructGlobal(frameId+1,blobId+1).surface.coeffs;
sigma1=zeros(3);
sigma1(1,1)=coeffs(1);sigma1(2,1)=coeffs(2);sigma1(3,1)=coeffs(3);
sigma1(2,2)=coeffs(4);sigma1(3,2)=coeffs(5);sigma1(3,3)=coeffs(6);
sigma1(1,2)=sigma1(2,1);sigma1(1,3)=sigma1(3,1);sigma1(2,3)=sigma1(3,2);

[U S V]=svd(sigma1,'econ');
S(3,3)=S(3,3)/sc;
sigma1=U*S*V';

blobStructGlobal(frameId+1,blobId+1).surface.coeffs(1:6)=[sigma1(1,1) sigma1(1,2) sigma1(1,3) sigma1(2,2) sigma1(2,3) sigma1(3,3)];


%delete current display of blob: we don't need to reload image
if(~isempty(handles.drawEllipseAxes123))
    
    for ii=1:length(handles.drawEllipseAxes123)
        if(handles.drawEllipseAxes123(ii)>0)%if object dissapears it returns 0
            delete(handles.drawEllipseAxes123(ii))
        end
    end
    handles.drawEllipseAxes123=[];
end

%draw new blob
ww=get(handles.sliderWindowRadius,'Value');
switch get(handles.popupmenuPyramidLevel,'Value')
    case 1
        ww=round([ww ww ww/handles.anisotropy]);
    case 2
        
        ww=round([ww ww 2.0*ww/handles.anisotropy]);
    case 3
        
    case 4
        
    case 5
        
    otherwise
end
cc=handles.centerTriViewAxes123;
handles.drawEllipseAxes123=drawEllipseOrthogonalPlanes(blobStructGlobal(frameId+1,blobId+1).surface.coeffs,...
    cc(1),cc(2),cc(3),cc(1)-ww(1),cc(2)-ww(2),cc(3)-ww(3),...
    handles.axes1XY,handles.axes2XZ,handles.axes3YZ,'g',get(handles.popupmenuPyramidLevel,'Value'));%redraw parent
drawnow();

%update the gui
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function sliderEditScaleZ_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sliderEditScaleZ (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

handles.sliderEditScaleZVal=get(hObject,'Value');%save value for later
%update the gui
guidata(hObject,handles);

% --- Executes on slider movement.
function sliderEditTranslateX_Callback(hObject, eventdata, handles)
% hObject    handle to sliderEditTranslateX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

global blobStructGlobal;

val=get(hObject,'Value');

if(handles.sliderEditTranslateXVal<val) %scale up
    sc=2.0;
else%scale down
    sc=-2.0;
end
handles.sliderEditTranslateXVal=val;%save value for later


%scale blob
blobId=get(handles.sliderBlob,'Value');
frameId=get(handles.sliderFrame,'Value');

if(isempty(blobStructGlobal(frameId+1,blobId+1).frame))
    set(handles.messageText,'String','Current object does not exist in this solution');
    return;
end

blobStructGlobal(frameId+1,blobId+1).surface.coeffs(7)=blobStructGlobal(frameId+1,blobId+1).surface.coeffs(7)+sc;
blobStructGlobal(frameId+1,blobId+1).center(1)=blobStructGlobal(frameId+1,blobId+1).center(1)+sc;


%delete current display of blob: we don't need to reload image
if(~isempty(handles.drawEllipseAxes123))
    
    for ii=1:length(handles.drawEllipseAxes123)
        if(handles.drawEllipseAxes123(ii)>0)%if object dissapears it returns 0
            delete(handles.drawEllipseAxes123(ii))
        end
    end
    handles.drawEllipseAxes123=[];
end

%draw new blob
ww=get(handles.sliderWindowRadius,'Value');
switch get(handles.popupmenuPyramidLevel,'Value')
    case 1
        ww=round([ww ww ww/handles.anisotropy]);
    case 2
        
        ww=round([ww ww 2.0*ww/handles.anisotropy]);
    case 3
        
    case 4
        
    case 5
        
    otherwise
end
cc=handles.centerTriViewAxes123;
handles.drawEllipseAxes123=drawEllipseOrthogonalPlanes(blobStructGlobal(frameId+1,blobId+1).surface.coeffs,...
    cc(1),cc(2),cc(3),cc(1)-ww(1),cc(2)-ww(2),cc(3)-ww(3),...
    handles.axes1XY,handles.axes2XZ,handles.axes3YZ,'g',get(handles.popupmenuPyramidLevel,'Value'));%redraw parent
drawnow();

%update the gui
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function sliderEditTranslateX_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sliderEditTranslateX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

handles.sliderEditTranslateXVal=get(hObject,'Value');%save value for later
%update the gui
guidata(hObject,handles);

% --- Executes on slider movement.
function sliderEditTranslateY_Callback(hObject, eventdata, handles)
% hObject    handle to sliderEditTranslateY (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

global blobStructGlobal;

val=get(hObject,'Value');

if(handles.sliderEditTranslateYVal<val) %scale up
    sc=2.0;
else%scale down
    sc=-2.0;
end
handles.sliderEditTranslateYVal=val;%save value for later


%scale blob
blobId=get(handles.sliderBlob,'Value');
frameId=get(handles.sliderFrame,'Value');

if(isempty(blobStructGlobal(frameId+1,blobId+1).frame))
    set(handles.messageText,'String','Current object does not exist in this solution');
    return;
end

blobStructGlobal(frameId+1,blobId+1).surface.coeffs(8)=blobStructGlobal(frameId+1,blobId+1).surface.coeffs(8)+sc;
blobStructGlobal(frameId+1,blobId+1).center(2)=blobStructGlobal(frameId+1,blobId+1).center(2)+sc;


%delete current display of blob: we don't need to reload image
if(~isempty(handles.drawEllipseAxes123))
    
    for ii=1:length(handles.drawEllipseAxes123)
        if(handles.drawEllipseAxes123(ii)>0)%if object dissapears it returns 0
            delete(handles.drawEllipseAxes123(ii))
        end
    end
    handles.drawEllipseAxes123=[];
end

%draw new blob
ww=get(handles.sliderWindowRadius,'Value');
switch get(handles.popupmenuPyramidLevel,'Value')
    case 1
        ww=round([ww ww ww/handles.anisotropy]);
    case 2
        
        ww=round([ww ww 2.0*ww/handles.anisotropy]);
    case 3
        
    case 4
        
    case 5
        
    otherwise
end
cc=handles.centerTriViewAxes123;
handles.drawEllipseAxes123=drawEllipseOrthogonalPlanes(blobStructGlobal(frameId+1,blobId+1).surface.coeffs,...
    cc(1),cc(2),cc(3),cc(1)-ww(1),cc(2)-ww(2),cc(3)-ww(3),...
    handles.axes1XY,handles.axes2XZ,handles.axes3YZ,'g',get(handles.popupmenuPyramidLevel,'Value'));%redraw parent
drawnow();


%update the gui
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function sliderEditTranslateY_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sliderEditTranslateY (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

handles.sliderEditTranslateYVal=get(hObject,'Value');%save value for later
%update the gui
guidata(hObject,handles);

% --- Executes on slider movement.
function sliderEditTranslateZ_Callback(hObject, eventdata, handles)
% hObject    handle to sliderEditTranslateZ (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

global blobStructGlobal;

val=get(hObject,'Value');

if(handles.sliderEditTranslateZVal<val) %scale up
    sc=2.0;
else%scale down
    sc=-2.0;
end
handles.sliderEditTranslateZVal=val;%save value for later


%scale blob
blobId=get(handles.sliderBlob,'Value');
frameId=get(handles.sliderFrame,'Value');

if(isempty(blobStructGlobal(frameId+1,blobId+1).frame))
    set(handles.messageText,'String','Current object does not exist in this solution');
    return;
end

blobStructGlobal(frameId+1,blobId+1).surface.coeffs(9)=blobStructGlobal(frameId+1,blobId+1).surface.coeffs(9)+sc;
blobStructGlobal(frameId+1,blobId+1).center(3)=blobStructGlobal(frameId+1,blobId+1).center(3)+sc;


%delete current display of blob: we don't need to reload image
if(~isempty(handles.drawEllipseAxes123))
    
    for ii=1:length(handles.drawEllipseAxes123)
        if(handles.drawEllipseAxes123(ii)>0)%if object dissapears it returns 0
            delete(handles.drawEllipseAxes123(ii))
        end
    end
    handles.drawEllipseAxes123=[];
end

%draw new blob
ww=get(handles.sliderWindowRadius,'Value');
switch get(handles.popupmenuPyramidLevel,'Value')
    case 1
        ww=round([ww ww ww/handles.anisotropy]);
    case 2
        
        ww=round([ww ww 2.0*ww/handles.anisotropy]);
    case 3
        
    case 4
        
    case 5
        
    otherwise
end
cc=handles.centerTriViewAxes123;
handles.drawEllipseAxes123=drawEllipseOrthogonalPlanes(blobStructGlobal(frameId+1,blobId+1).surface.coeffs,...
    cc(1),cc(2),cc(3),cc(1)-ww(1),cc(2)-ww(2),cc(3)-ww(3),...
    handles.axes1XY,handles.axes2XZ,handles.axes3YZ,'g',get(handles.popupmenuPyramidLevel,'Value'));%redraw parent
drawnow();

%update the gui
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function sliderEditTranslateZ_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sliderEditTranslateZ (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

handles.sliderEditTranslateZVal=get(hObject,'Value');%save value for later
%update the gui
guidata(hObject,handles);

% --- Executes on slider movement.
function slider17_Callback(hObject, eventdata, handles)
% hObject    handle to slider17 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


%update the gui
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function slider17_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider17 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


%update the gui
guidata(hObject,handles);

% --- Executes on slider movement.
function sliderEditTheta_Callback(hObject, eventdata, handles)
% hObject    handle to sliderEditTheta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

global blobStructGlobal;

val=get(hObject,'Value');

if(handles.sliderEditThetaVal<val) %scale up
    sc=15;
else%scale down
    sc=-15;
end
handles.sliderEditThetaVal=val;%save value for later


%scale blob
blobId=get(handles.sliderBlob,'Value');
frameId=get(handles.sliderFrame,'Value');

if(isempty(blobStructGlobal(frameId+1,blobId+1).frame))
    set(handles.messageText,'String','Current object does not exist in this solution');
    return;
end

coeffs=blobStructGlobal(frameId+1,blobId+1).surface.coeffs;
sigma1=zeros(3);
sigma1(1,1)=coeffs(1);sigma1(2,1)=coeffs(2);sigma1(3,1)=coeffs(3);
sigma1(2,2)=coeffs(4);sigma1(3,2)=coeffs(5);sigma1(3,3)=coeffs(6);
sigma1(1,2)=sigma1(2,1);sigma1(1,3)=sigma1(3,1);sigma1(2,3)=sigma1(3,2);

[V D]=eig(sigma1);
R=[cosd(sc) -sind(sc) 0;...
   sind(sc) cosd(sc)  0;
   0 0 1];

%reorganize V to make sure we rotate around XY
[aux idx]=sort(abs(V(3,:)));
V=V(:,idx);
aux=diag(D);
D=diag(aux(idx));
V=V*R;

sigma1=V*D*V';

blobStructGlobal(frameId+1,blobId+1).surface.coeffs(1:6)=[sigma1(1,1) sigma1(1,2) sigma1(1,3) sigma1(2,2) sigma1(2,3) sigma1(3,3)];


%delete current display of blob: we don't need to reload image
if(~isempty(handles.drawEllipseAxes123))
    
    for ii=1:length(handles.drawEllipseAxes123)
        if(handles.drawEllipseAxes123(ii)>0)%if object dissapears it returns 0
        delete(handles.drawEllipseAxes123(ii))
        end
    end
    handles.drawEllipseAxes123=[];
end

%draw new blob
ww=get(handles.sliderWindowRadius,'Value');
switch get(handles.popupmenuPyramidLevel,'Value')
    case 1
        ww=round([ww ww ww/handles.anisotropy]);
    case 2
        
        ww=round([ww ww 2.0*ww/handles.anisotropy]);
    case 3
        
    case 4
        
    case 5
        
    otherwise
end
cc=handles.centerTriViewAxes123;
handles.drawEllipseAxes123=drawEllipseOrthogonalPlanes(blobStructGlobal(frameId+1,blobId+1).surface.coeffs,...
    cc(1),cc(2),cc(3),cc(1)-ww(1),cc(2)-ww(2),cc(3)-ww(3),...
    handles.axes1XY,handles.axes2XZ,handles.axes3YZ,'g',get(handles.popupmenuPyramidLevel,'Value'));%redraw parent
drawnow();

%update the gui
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function sliderEditTheta_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sliderEditTheta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

handles.sliderEditThetaVal=get(hObject,'Value');%save value for later
%update the gui
guidata(hObject,handles);

% --- Executes on slider movement.
function sliderEditPhi_Callback(hObject, eventdata, handles)
% hObject    handle to sliderEditPhi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


%update the gui
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function sliderEditPhi_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sliderEditPhi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

handles.sliderEditPhiVal=get(hObject,'Value');%save value for later
%update the gui
guidata(hObject,handles);

% --- Executes on button press in pushbuttonSaveSolution.
function pushbuttonSaveSolution_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonSaveSolution (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

saveFile=['./saveTemp/blobStruct_stackMCMC_' datestr(now,30) '.mat'];
save(saveFile,'-struct','handles','blobStruct','stackFilename');
set(handles.messageText,'String',['Solution saved temporarily at ' saveFile]);



%update the gui
guidata(hObject,handles);

%------------------------------------------------------------------------
%If we have loaded a list with some sort of score for this frame it
%displays the object in order
% --- Executes on button press in pushbuttonNextSorted.
function pushbuttonNextSorted_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonNextSorted (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles=pushbuttonSorted(hObject,eventdata,handles,true);
guidata(hObject,handles);


%dso I can call it from other functions and handles information is
%preserved
function handles=pushbuttonSorted(hObject,eventdata,handles,isForward)

global blobStructGlobal;

frameId=get(handles.sliderFrame,'Value');

if(handles.nextObjectSorted(2)~=frameId)%we need to recalculate sorted list
    set(handles.messageText,'String','Recalculating sorted list for current frame');
    %------------------change this depending on what criteria you want to
    %use---------------------------------------
    %[handles.nextObjectSortedList, handles.nextObjectSortedListVal] =sortFrameElementsBySize(blobStructGlobal{frameId+1},[1 1 handles.anisotropy]);
    %[handles.nextObjectSortedList, handles.nextObjectSortedListVal] =sortFrameElementsByNumberOfSupervoxels(blobStructGlobal{frameId+1});
    %[handles.nextObjectSortedList, handles.nextObjectSortedListVal] = sortFrameElementsByScore(blobStructGlobal{frameId+1});
    %[handles.nextObjectSortedList, handles.nextObjectSortedListVal] = sortFrameElementsByBetaPrior(blobStructGlobal{frameId+1});%used to store probBAckground
    %[handles.nextObjectSortedList, handles.nextObjectSortedListVal] = sortFrameByLineageLength(frameId);
    %[handles.nextObjectSortedList, handles.nextObjectSortedListVal] =sortFrameElementsByEccentricity(blobStructGlobal{frameId+1},[1 1 handles.anisotropy]);
    %main one for cell division
    [handles.nextObjectSortedList, handles.nextObjectSortedListVal] = sortFrameElementsByProbCellDivision(blobStructGlobal{frameId+1},blobStructGlobal{frameId+2});
    %[handles.nextObjectSortedList, handles.nextObjectSortedListVal] =sortFrameElementsByCellDivision(blobStructGlobal{frameId+1},blobStructGlobal{frameId+2});
    %[handles.nextObjectSortedList, handles.nextObjectSortedListVal] =sortFrameElementsByTrackBirth(blobStructGlobal{frameId+1});
    %[handles.nextObjectSortedList, handles.nextObjectSortedListVal] =sortFrameElementsByDisplacementWithRespectToParent(blobStructGlobal{frameId+1},blobStructGlobal{frameId+2},handles.anisotropy);
    %[handles.nextObjectSortedList, handles.nextObjectSortedListVal] = sortFrameElementsByZvalue(blobStructGlobal{frameId+1});
    %[handles.nextObjectSortedList, handles.nextObjectSortedListVal] =sortFrameElementsBySizeChangeWithRespectToParent(blobStructGlobal{frameId+1},blobStructGlobal{frameId+2},handles.anisotropy);
    %[handles.nextObjectSortedList, handles.nextObjectSortedListVal] =sortFrameElementsByCellDeath(blobStructGlobal{frameId+1});
    %[handles.nextObjectSortedList, handles.nextObjectSortedListVal] =sortFrameElementsByLineageBirth(blobStructGlobal{frameId+1});
    %[handles.nextObjectSortedList, handles.nextObjectSortedListVal] =[1:size(blobStructGlobal,2)];%no sroting. just in order
    %use precomputed order
    %load('/Users/amatf/TrackingNuclei/tmp/GMEMtracking3D_1316463792_GPU_thr002_splitAdaBoostm70maskSph/accuracyStats/cellDivisionDetection-adaBoost/sortOrderCellDivision.mat','sortOrder');
    %[handles.nextObjectSortedList, handles.nextObjectSortedListVal] =sortOrder;
    %--------------------------------------------
    handles.nextObjectSorted(2)=frameId;
    handles.nextObjectSorted(1)=0;
    set(handles.messageText,'String','Done recalculating sorted list for current frame');
end
    
if(isForward == true )
    handles.nextObjectSorted(1) = handles.nextObjectSorted(1)+1;
else
    handles.nextObjectSorted(1) = handles.nextObjectSorted(1)-1;
end

if(handles.nextObjectSorted(1)<=0)
    handles.nextObjectSorted(1)=0;
    
    set(handles.messageText,'String','No more previous cell divisions in this frame');
               
    guidata(hObject,handles);
    return;
elseif(handles.nextObjectSorted(1)>length(handles.nextObjectSortedList) || handles.nextObjectSortedListVal(handles.nextObjectSorted(1)) < 0)
    set(handles.messageText,'String','No more cell divisions in this frame');
    
    handles.nextObjectSorted(1) = handles.nextObjectSorted(1)-1;
    %handles.nextObjectSorted(1)=1;
        
    guidata(hObject,handles);
    return;
    
end

blobId = handles.nextObjectSortedList(handles.nextObjectSorted(1))-1;

disp(['Sorted val = ' num2str(handles.nextObjectSortedListVal(handles.nextObjectSorted(1))) ])

%disp 'WARNING:patch to annotate'
%blobId=get(handles.sliderBlob,'Value')+1;
%handles.classifierResultsSingleFrame(blobId+1)
set(handles.sliderBlob,'Value',blobId);

%--------------
% % disp(['ColorIdx=' num2str(blobStructGlobal(frameId+1,blobId+1).colorIdx)])
% % if(blobStructGlobal(frameId+1,blobId+1).colorIdx<5)
% %     error 'No more divisions!!!'
% % end
%----------------

guidata(hObject,handles);
sliderBlob_Callback(hObject, eventdata, handles);

%guidata(hObject,handles);%slider callback already saves handles
%---------------------------------------------------------------------

% --- Executes on button press in pushbuttonEditDeleteTrack.
function pushbuttonEditDeleteTrack_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonEditDeleteTrack (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global blobStructGlobal;

blobId=get(handles.sliderBlob,'Value');
frameId=get(handles.sliderFrame,'Value');

if(isempty(blobStructGlobal(frameId+1,blobId+1).frame) || blobStructGlobal(frameId+1,blobId+1).center(1)<-1e31)
    return;
end


ff=frameId+1;
bb=blobId+1;

%---------------------------------------
%{
%A-go to the top of the tree
while(blobStructGlobal(ff,bb).solutions(1).parentIdx(1)<4e9)
    ffAux=blobStructGlobal(ff,bb).solutions(1).parentIdx(1)+1;
    bb=blobStructGlobal(ff,bb).solutions(1).parentIdx(2)+1;
    ff=ffAux;
end
%}

%B-start at the selected nuclei and only go down the tree
if(blobStructGlobal(ff,bb).solutions(1).parentIdx(1)<4e9)
    ffAux=blobStructGlobal(ff,bb).solutions(1).parentIdx(1)+1;
    bbAux=blobStructGlobal(ff,bb).solutions(1).parentIdx(2)+1;
    
    ch=blobStructGlobal(ffAux,bbAux).solutions(1).childrenIdx;
    pp=find(ch==bb-1);
    if(length(pp)~=1 || pp<2)
        error 'Tree is wrong'
    else        
        blobStructGlobal(ffAux,bbAux).solutions(1).childrenIdx([pp pp-1])=[];
    end
end
%------------------------------------------------


%traverse the tree down
queue=[ff bb];
while(~isempty(queue))
    ff=queue(1,1);
    bb=queue(1,2);
    queue(1,:)=[];
    
    ch=blobStructGlobal(ff,bb).solutions(1).childrenIdx;    
    center=blobStructGlobal(ff,bb).surface.coeffs(7:9);
    
    if(center(1)<-1e31)%dead cell: it can not have children
       continue; 
    end
    
    for ii=1:2:length(ch)
        if(length(ch)>2)%cell division
            queue=[queue;ch(ii)+1 ch(ii+1)+1];
        else
            queue=[queue;ch(ii)+1 ch(ii+1)+1];
        end
    end
    %kill cell
    blobStructGlobal(ff,bb).center=repmat(-1e32,size(blobStructGlobal(ff,bb).center));
    blobStructGlobal(ff,bb).surface.coeffs(7:9)=repmat(-1e32,size(blobStructGlobal(ff,bb).center));
    blobStructGlobal(ff,bb).solutions(1).parentIdx=[uint32(2^32-1) 0];
    blobStructGlobal(ff,bb).solutions(1).childrenIdx=[];
    blobStructGlobal(ff,bb).solutions(1).label=0;
end

set(handles.messageText,'String','Track deleted');
%update the gui
guidata(hObject,handles);


% --- Executes on button press in checkboxLineageDisplay.
function checkboxLineageDisplay_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxLineageDisplay (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkboxLineageDisplay



function editZanisotropy_Callback(hObject, eventdata, handles)
% hObject    handle to editZanisotropy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editZanisotropy as text
%        str2double(get(hObject,'String')) returns contents of editZanisotropy as a double

handles.anisotropy=str2double(get(hObject,'String'));
%update the gui
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function editZanisotropy_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editZanisotropy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

handles.anisotropy=str2double(get(hObject,'String'));
%update the gui
guidata(hObject,handles);


% --- Executes on selection change in popupmenuPyramidLevel.
function popupmenuPyramidLevel_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenuPyramidLevel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenuPyramidLevel contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenuPyramidLevel

frameId = get(handles.sliderFrame,'Value');

handles=cacheStackImages(handles,frameId+1,hObject);

pathname = handles.pathLogFile;
pathname = [pathname 'XML_finalResult_lht' filesep];
handles.svFilename = cacheStackSupervoxels(handles, frameId+1, pathname);



cacheTGMMxmlSolution(pathname, frameId,handles.frameIni, 0 );%TODO: fix nearest neighbors value from GUI

set(handles.messageText,'String',['Please point to a new blob or frame to upload triview scale']);

%update the gui
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function popupmenuPyramidLevel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenuPyramidLevel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

set(hObject,'Value',1);%default value
%radius: we load +-stackMaxSizeCache
%update the gui
guidata(hObject,handles);


% --- Executes on selection change in popupmenuSpecializedFunctions.
function popupmenuSpecializedFunctions_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenuSpecializedFunctions (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenuSpecializedFunctions contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenuSpecializedFunctions

global stackGlobal;
global blobStructGlobal;

blobId=get(handles.sliderBlob,'Value');
frameId=get(handles.sliderFrame,'Value');


numNN=15;%number of nearest neighbors to display
deltaT=4;%+- time to display lineages
switch(get(hObject,'Value'))
    
    case 1%Neigboring lineage visulization
        rectCoord=[1e32 0 1e32 0 1e32 0];%xMin xMax yMin yMax zMin zMax
        knn=blobStructGlobal(frameId+1,blobId+1).neigh;
        for kk=1:size(knn,1)
           if(knn(kk,1)~=frameId) continue;end;
           cc=blobStructGlobal(frameId+1,knn(kk,2)+1).center;
           if(cc(1)<-1e31) continue;end;
           
           rectCoord(1)=min(cc(1),rectCoord(1));
           rectCoord(2)=max(cc(1),rectCoord(2));
           rectCoord(3)=min(cc(2),rectCoord(3));
           rectCoord(4)=max(cc(2),rectCoord(4));
           %rectCoord(5)=min(cc(3),rectCoord(5));
           %rectCoord(6)=max(cc(3),rectCoord(6));
        end
        imSize=size(stackGlobal{frameId+1-handles.stackIniCache+1});
        rectCoord(1:2:end)=rectCoord(1:2:end)-30;
        rectCoord(2:2:end)=rectCoord(2:2:end)+30;%add some margin
        rectCoord([5 6])=[blobStructGlobal(frameId+1,blobId+1).center(3)-10 blobStructGlobal(frameId+1,blobId+1).center(3)+10];%define Z coordinate
               
        
        
        switch get(handles.popupmenuPyramidLevel,'Value')
            case 1
                
            case 2
                rectCoord(1:4)=(rectCoord(1:4)+1)/2.0;%+1 is vecause of Matlab indexing in pyramid. 1:2:end->pixel3 goes to pixel2=0.5*(3+1) in next pyramid level                
            case 3
                rectCoord(1:4)=(rectCoord(1:4)+1)/4.0;%in the first two levels only x,y are downsampled                
            case 4
                rectCoord(1:4)=(rectCoord(1:4)+1)/8.0;
                rectCoord([5 6])=(rectCoord([5 6])+1)/2.0;                
            case 5
                rectCoord(1:4)=(rectCoord(1:4)+1)/16.0;
                rectCoord([5 6])=(rectCoord([5 6])+1)/4.0;                
            otherwise
        end
        %set boundaries
        rectCoord=round(rectCoord);
        rectCoord(1)=max(1,rectCoord(1));
        rectCoord(2)=min(imSize(1),rectCoord(2));
        rectCoord(3)=max(1,rectCoord(3));
        rectCoord(4)=min(imSize(2),rectCoord(4));
        rectCoord(5)=max(1,rectCoord(5));
        rectCoord(6)=min(imSize(2),rectCoord(6));
        
        %generate image structure to display local lineage tree
        imStruct(1).rectCoord=rectCoord;
        imStruct(2*deltaT+1).rectCoord=imStruct(1).rectCoord;%preallocate memory
        minIm=1e32;maxIm=0;%to scale images
        for kk=1:2*deltaT+1
            imStruct(kk).rectCoord=imStruct(1).rectCoord;
            imStruct(kk).stack=stackGlobal{frameId+1+kk-1-deltaT-handles.stackIniCache+1}(rectCoord(1):rectCoord(2),rectCoord(3):rectCoord(4),rectCoord(5):rectCoord(6));
            minIm=min(minIm,min(imStruct(kk).stack(:)));
            maxIm=max(maxIm,max(imStruct(kk).stack(:)));
        end
        minIm=single(minIm);maxIm=single(maxIm);
        for kk=1:2*deltaT+1
           stackAux=single(imStruct(kk).stack);
           imStruct(kk).stack=uint8(255.0*(stackAux-minIm)/(maxIm-minIm));
        end
        
        localLineageDisplay(blobStructGlobal(frameId+1-deltaT:frameId+1+deltaT,:),frameId+1,blobId+1,numNN,deltaT,handles.anisotropy,imStruct,frameId-deltaT,get(handles.popupmenuPyramidLevel,'Value'));%display 10 nearest neighbors
    otherwise
        a=1;%do nothing
end

%update the gui
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function popupmenuSpecializedFunctions_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenuSpecializedFunctions (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbuttonCorrectCellDivision.
function pushbuttonCorrectCellDivision_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonCorrectCellDivision (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%saves 3 annotations at once: mother + daughter1+daughter2
handles=saveAnnotationInMemorySupervoxels(handles,hObject,3);
guidata(hObject,handles);


% --- Executes on button press in pushbuttonShouldBeCellDivision.
function pushbuttonShouldBeCellDivision_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonShouldBeCellDivision (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%saves 3 annotations at once: mother + daughter1+daughter2
handles=saveAnnotationInMemorySupervoxels(handles,hObject,4);
guidata(hObject,handles);


% --- Executes on button press in pushbuttonWrongCellDivision.
function pushbuttonWrongCellDivision_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonWrongCellDivision (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%saves 3 annotations at once: mother + daughter1+daughter2
handles=saveAnnotationInMemorySupervoxels(handles,hObject,5);
guidata(hObject,handles);



function editMaxCacheSize_Callback(hObject, eventdata, handles)
% hObject    handle to editMaxCacheSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editMaxCacheSize as text
%        str2double(get(hObject,'String')) returns contents of editMaxCacheSize as a double

handles.stackMaxSizeCache = str2double(get(hObject,'String'));
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function editMaxCacheSize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editMaxCacheSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkboxShowBlobSupervoxels.
function checkboxShowBlobSupervoxels_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxShowBlobSupervoxels (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkboxShowBlobSupervoxels

if( get(hObject,'Value') == 0 )
    %make everything invisible
     if( isempty(handles.drawSupervoxelAxes123) == false)
         for kk = 1:numel(handles.drawSupervoxelAxes123)
            if( handles.drawSupervoxelAxes123(kk) >= 0 )
               set( handles.drawSupervoxelAxes123(kk), 'Visible','off');
            end
         end
     end
    guidata(hObject,handles);
else %make supervoxel visible
    if( isempty(handles.drawSupervoxelAxes123) )
        %redraw everything
        sliderBlob_Callback(hObject, eventdata, handles);
    else
        %just make everything visible
        for kk = 1:numel(handles.drawSupervoxelAxes123)
            if( handles.drawSupervoxelAxes123(kk) >= 0 )
               set( handles.drawSupervoxelAxes123(kk), 'Visible','on');
            end
        end
        guidata(hObject,handles); 
    end
end

% --- Executes on button press in checkboxShowAllSupervoxels.
function checkboxShowAllSupervoxels_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxShowAllSupervoxels (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkboxShowAllSupervoxels


% --- Executes on button press in checkbox4.
function checkbox4_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox4


% --- Executes on button press in pushbuttonTriviewXup.
function pushbuttonTriviewXup_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonTriviewXup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


handles.offsetCenterPlane123(1) = handles.offsetCenterPlane123(1) + 1;
guidata(hObject,handles);
sliderBlob_Callback(hObject, eventdata, handles);

% --- Executes on button press in pushbuttonTriviewXdown.
function pushbuttonTriviewXdown_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonTriviewXdown (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.offsetCenterPlane123(1) = handles.offsetCenterPlane123(1) - 1;
guidata(hObject,handles);
sliderBlob_Callback(hObject, eventdata, handles);


% --- Executes on button press in pushbuttonTriviewYup.
function pushbuttonTriviewYup_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonTriviewYup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.offsetCenterPlane123(2) = handles.offsetCenterPlane123(2) + 1;
guidata(hObject,handles);
sliderBlob_Callback(hObject, eventdata, handles);


% --- Executes on button press in pushbuttonTriviewYdown.
function pushbuttonTriviewYdown_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonTriviewYdown (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.offsetCenterPlane123(2) = handles.offsetCenterPlane123(2) - 1;
guidata(hObject,handles);
sliderBlob_Callback(hObject, eventdata, handles);


% --- Executes on button press in pushbuttonTriviewZup.
function pushbuttonTriviewZup_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonTriviewZup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.offsetCenterPlane123(3) = handles.offsetCenterPlane123(3) + 1;
guidata(hObject,handles);
sliderBlob_Callback(hObject, eventdata, handles);


% --- Executes on button press in pushbuttonZdown.
function pushbuttonZdown_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonZdown (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.offsetCenterPlane123(3) = handles.offsetCenterPlane123(3) - 1;
guidata(hObject,handles);
sliderBlob_Callback(hObject, eventdata, handles);


% --- Executes on button press in pushbuttonTriviewRecenter.
function pushbuttonTriviewRecenter_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonTriviewRecenter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.offsetCenterPlane123 = [0 0 0];
guidata(hObject,handles);
sliderBlob_Callback(hObject, eventdata, handles);


% --- Executes on button press in pushbuttonTriviewXup456.
function pushbuttonTriviewXup456_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonTriviewXup456 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.offsetCenterPlane456(1) = handles.offsetCenterPlane456(1) + 1;
guidata(hObject,handles);
sliderBlob_Callback(hObject, eventdata, handles);

% --- Executes on button press in pushbuttonTriviewXdown456.
function pushbuttonTriviewXdown456_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonTriviewXdown456 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.offsetCenterPlane456(1) = handles.offsetCenterPlane456(1) - 1;
guidata(hObject,handles);
sliderBlob_Callback(hObject, eventdata, handles);

% --- Executes on button press in pushbuttonTriviewYup456.
function pushbuttonTriviewYup456_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonTriviewYup456 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.offsetCenterPlane456(2) = handles.offsetCenterPlane456(2) + 1;
guidata(hObject,handles);
sliderBlob_Callback(hObject, eventdata, handles);

% --- Executes on button press in pushbuttonTriviewYdown456.
function pushbuttonTriviewYdown456_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonTriviewYdown456 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.offsetCenterPlane456(2) = handles.offsetCenterPlane456(2) - 1;
guidata(hObject,handles);
sliderBlob_Callback(hObject, eventdata, handles);

% --- Executes on button press in pushbuttonTriviewZup456.
function pushbuttonTriviewZup456_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonTriviewZup456 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.offsetCenterPlane456(3) = handles.offsetCenterPlane456(3) + 1;
guidata(hObject,handles);
sliderBlob_Callback(hObject, eventdata, handles);

% --- Executes on button press in pushbuttonTriviewZdown456.
function pushbuttonTriviewZdown456_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonTriviewZdown456 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.offsetCenterPlane456(3) = handles.offsetCenterPlane456(3) - 1;
guidata(hObject,handles);
sliderBlob_Callback(hObject, eventdata, handles);

% --- Executes on button press in pushbuttonTriviewRecenter456.
function pushbuttonTriviewRecenter456_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonTriviewRecenter456 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.offsetCenterPlane456 = [0 0 0];
guidata(hObject,handles);
sliderBlob_Callback(hObject, eventdata, handles);



function edit5numNN_Callback(hObject, eventdata, handles)
% hObject    handle to edit5numNN (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit5numNN as text
%        str2double(get(hObject,'String')) returns contents of edit5numNN as a double


% --- Executes during object creation, after setting all properties.
function edit5numNN_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit5numNN (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbuttonPrevSorted.
function pushbuttonPrevSorted_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonPrevSorted (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


handles=pushbuttonSorted(hObject,eventdata,handles,false);
guidata(hObject,handles);


% --- Executes on button press in pushbuttonGenFeatures.
function pushbuttonGenFeatures_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonGenFeatures (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%collect features from the tGMM run by matching XYZ coordinates
set(handles.messageText,'String',['Collecting features from annotations and TGMM run...']);
%read annotations
pathAnn = [handles.pathLogFile 'annForCellDivDiscrWithTempWin'];
addpath([fileparts( mfilename('fullpath') ) filesep '..' filesep ]);
[obj, yClass] = parseXmlAnnotationsFolder(pathAnn);
addpath([fileparts( mfilename('fullpath') ) filesep '..' filesep ]);

%parse xyz and TM for each annotation
posTM = findstr(handles.imgFilePattern, '?');
pAux = find(diff(posTM) > 1);
if( isempty(pAux) == false )
   posTM(pAux(1)+1:end) = []; 
end

annN = length(obj);
xyzt = zeros(annN, 4);

for ii = 1:annN
    xyzt(ii,:) = [obj(ii).m, str2double(obj(ii).imFilename(posTM))];
end

%read features for each annotation
addpath([fileparts( mfilename('fullpath') ) filesep '..' filesep ]);
xTrainAll = [];
yTrainAll = [];
for TM = unique(xyzt(:,4))'    
    pos = find(xyzt(:,4) == TM);    
    xyzAnn = xyzt(pos,1:3);
    yLabel = yClass(pos);
    %read xyz coordinate txt file
    xyzF = load([handles.pathLogFile 'CDTWfeatures' filesep 'CDWTfeatures_TM' num2str(TM,'%.5d') '.txt']);
    
    %read binary file with features
    [xTrain, yTrain] = readTrainingDataBinary([handles.pathLogFile 'CDTWfeatures' filesep 'CDWTfeatures_TM' num2str(TM,'%.5d') '.bin']);
    
    %find the match
    idx = knnsearch(xyzF, xyzAnn, 'k',1);
    
    xTrainAll = [xTrainAll; xTrain(idx,:)];
    yTrainAll = [yTrainAll; yLabel];             
end

%write features
outFilename = [handles.pathLogFile 'annForCellDivDiscrWithTempWin' filesep 'trainFeaturesCDTW.bin'];
writeTrainingDataBinary(xTrainAll, yTrainAll, outFilename);

rmpath([fileparts( mfilename('fullpath') ) filesep '..' filesep ]);

set(handles.messageText,'String',['All features saved successfully at ' outFilename]);


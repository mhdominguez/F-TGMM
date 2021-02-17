function varargout = trainClassifier(varargin)
% TRAINCLASSIFIER MATLAB code for trainClassifier.fig
%      TRAINCLASSIFIER, by itself, creates a new TRAINCLASSIFIER or raises the existing
%      singleton*.
%
%      H = TRAINCLASSIFIER returns the handle to a new TRAINCLASSIFIER or the handle to
%      the existing singleton*.
%
%      TRAINCLASSIFIER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in TRAINCLASSIFIER.M with the given input arguments.
%
%      TRAINCLASSIFIER('Property','Value',...) creates a new TRAINCLASSIFIER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before trainClassifier_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to trainClassifier_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help trainClassifier

% Last Modified by GUIDE v2.5 06-Nov-2014 16:44:28

% Begin initialization code - DO NOT EDIT
gui_Singleton = 0;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @trainClassifier_OpeningFcn, ...
                   'gui_OutputFcn',  @trainClassifier_OutputFcn, ...
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


% --- Executes just before trainClassifier is made visible.
function trainClassifier_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to trainClassifier (see VARARGIN)

global xTrainGlobal;
global yTrainGlobal;


handles.midPlaneFeatureIndex = 1027 + 1;%you have to make sure this agrees with the C++ code generation

% Choose default command line output for trainClassifier
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes trainClassifier wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = trainClassifier_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbuttonTrainClassifier.
function pushbuttonTrainClassifier_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonTrainClassifier (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global xTrainGlobal;
global yTrainGlobal;

thr = str2double( get(handles.edit4ThrCellDivMidPlane,'string')) + 0.5;%to give some room
pos = find( xTrainGlobal(:,handles.midPlaneFeatureIndex) < thr );




kCV = 10;
minLeaf = str2double( get(handles.editMinElemPerLeaf,'string'));
numWeakLearners = str2double( get(handles.editNumWeakLearners,'string'));
learnRate = str2double( get(handles.editLearnRate,'string'));


set(handles.text2,'String','Calculating boosting tree classifier...it can take a while');
set(handles.text2,'BackgroundColor',[0 1 1])
drawnow();

disp('=========TODO: add symmetry as a GUI parameter====================');
symmetry = 1;
[rusTree, thrVec, prec, rec] = trainClassifierRusboost( xTrainGlobal(pos,:), yTrainGlobal(pos), symmetry, numWeakLearners, minLeaf, learnRate, kCV, handles.axes2);

set(handles.text2,'String','Classifier successfully calculated');
set(handles.text2,'BackgroundColor',[0 1 0]);



plot(handles.axes1, thrVec,rec,'k')
xlabel(handles.axes1,'Threshold classifier');
ylabel(handles.axes1,'Recall');
grid(handles.axes1,'on');

hold(handles.axes1,'on');
plot(handles.axes1, thrVec,prec,'g')
xlabel(handles.axes1,'Threshold classifier');
ylabel(handles.axes1,'Precision');
grid(handles.axes1,'on');
hold(handles.axes1,'off');

legend(handles.axes1, 'Recall', 'Precision','Location','Best');




%save changes
guidata(hObject,handles);


function editNumWeakLearners_Callback(hObject, eventdata, handles)
% hObject    handle to editNumWeakLearners (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editNumWeakLearners as text
%        str2double(get(hObject,'String')) returns contents of editNumWeakLearners as a double


% --- Executes during object creation, after setting all properties.
function editNumWeakLearners_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editNumWeakLearners (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editMinElemPerLeaf_Callback(hObject, eventdata, handles)
% hObject    handle to editMinElemPerLeaf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editMinElemPerLeaf as text
%        str2double(get(hObject,'String')) returns contents of editMinElemPerLeaf as a double


% --- Executes during object creation, after setting all properties.
function editMinElemPerLeaf_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editMinElemPerLeaf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editLearnRate_Callback(hObject, eventdata, handles)
% hObject    handle to editLearnRate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editLearnRate as text
%        str2double(get(hObject,'String')) returns contents of editLearnRate as a double


% --- Executes during object creation, after setting all properties.
function editLearnRate_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editLearnRate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbuttonTrainFinalClassifier.
function pushbuttonTrainFinalClassifier_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonTrainFinalClassifier (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global xTrainGlobal;
global yTrainGlobal;

thr = str2double( get(handles.edit4ThrCellDivMidPlane,'string')) + 0.5;%to give some room
pos = find( xTrainGlobal(:,handles.midPlaneFeatureIndex) < thr );

minLeaf = str2double( get(handles.editMinElemPerLeaf,'string'));
numWeakLearners = str2double( get(handles.editNumWeakLearners,'string'));
learnRate = str2double( get(handles.editLearnRate,'string'));



set(handles.text2,'String','Calculating final boosting tree classifier...it can take a while');
set(handles.text2,'BackgroundColor',[0 1 1])
drawnow();

t = ClassificationTree.template('minleaf',minLeaf);
rusTree = fitensemble(xTrainGlobal(pos,:), yTrainGlobal(pos),'RUSBoost',numWeakLearners,t,'LearnRate',learnRate,'nprint',20);

%write classifier
filenameClassifier = [handles.pathFeatureFiles 'classifierCellDivisionWithTemporalWindow.txt'];
parseMatlabFitEnsembleToCppGentleBoostFormat(rusTree, filenameClassifier);

%write classifier as a Matlab file
filenameClassifier = [handles.pathFeatureFiles 'classifierCellDivisionWithTemporalWindow.mat'];
save(filenameClassifier, 'rusTree');

set(handles.text2,'String',['Classifier successfully calculated. Saved to file ' filenameClassifier '. Copy it to TGMM bin folder if you want to use it for tracking.']);
set(handles.text2,'BackgroundColor',[0 1 0]);

%save changes
guidata(hObject,handles);


function edit4ThrCellDivMidPlane_Callback(hObject, eventdata, handles)
% hObject    handle to edit4ThrCellDivMidPlane (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit4ThrCellDivMidPlane as text
%        str2double(get(hObject,'String')) returns contents of edit4ThrCellDivMidPlane as a double


% --- Executes during object creation, after setting all properties.
function edit4ThrCellDivMidPlane_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit4ThrCellDivMidPlane (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbuttonAdjustMidPlaneThr.
function pushbuttonAdjustMidPlaneThr_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonAdjustMidPlaneThr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global xTrainGlobal;
global yTrainGlobal;

%thr = str2double( get(handles.edit4ThrCellDivMidPlane,'string'));
xTrain = xTrainGlobal(:,handles.midPlaneFeatureIndex);

[thrVec,prec,rec] = fitMidplaneThr(xTrain, yTrainGlobal,5);
plot(handles.axes1, thrVec,mean(rec,2))
xlabel(handles.axes1,'Threshold Cell Div. Plane value');
ylabel(handles.axes1,'Recall');
grid(handles.axes1,'on');
set(handles.text2,'String','Adjust the thr. Cell Div Plane value base on recall. This is a first filter so threshold should be set to very high recall to not lose cell division.');
set(handles.text2,'BackgroundColor',[1 1 0])

%just in case user forgets
thr = thrVec(mean(rec,2) > 0.99);
set(handles.edit4ThrCellDivMidPlane,'string',num2str(thr(1)));

%save changes
guidata(hObject,handles);

% --- Executes on button press in pushbuttonLoadFeatures.
function pushbuttonLoadFeatures_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonLoadFeatures (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%obtain path to stack
[filename,pathname] = uigetfile('*.bin', 'Select path binary file with features');
if(filename==0)
    return;
end

handles.pathFeatureFiles = pathname;

D = dir([pathname '*.bin']);

global xTrainGlobal;
global yTrainGlobal;
%read training data
xTrainGlobal = [];
yTrainGlobal = [];
for ii = 1:length(D)
    binaryFeatureFiles = [pathname D(ii).name];
  [xTrainAux, yTrainAux] = readTrainingDataBinary(binaryFeatureFiles);  
  xTrainGlobal = [ xTrainGlobal; xTrainAux];
  yTrainGlobal = [ yTrainGlobal; yTrainAux];
end

clear xTrainAux yTrainAux

set(handles.text2,'String',[num2str(length(yTrainGlobal)) ' samples loaded successfully']);
set(handles.text2,'BackgroundColor',[0 1 0])

%save changes
guidata(hObject,handles);

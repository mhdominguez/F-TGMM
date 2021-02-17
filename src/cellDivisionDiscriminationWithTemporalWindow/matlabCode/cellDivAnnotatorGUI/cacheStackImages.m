function handles=cacheStackImages(handles,centerImgIdx,hObject)
global stackGlobal;

set(handles.messageText,'String',['Loading frames...']);
guidata(hObject,handles);
drawnow();

maxSize=handles.stackMaxSizeCache;
numFrames=length(handles.stackFilename);
imgFilename=handles.stackFilename;




%check if we need to read tif or jp2
readType = 0;%0->tif;1->jp2
fileExtension = '.tif';
for kk = 1:length(imgFilename)
   if(isempty(imgFilename{kk}) == false)
       if ( strcmp( imgFilename{kk}(end-2:end), 'jp2' ) == true)
           readType = 1;
           fileExtension = '.jp2';
       elseif ( strcmp( imgFilename{kk}(end-2:end), 'klb' ) == true)
           readType = 2;
           fileExtension = '.klb';
       end
       
       break;
   end
end


%generate filenames
%outside because parfor does not like handles or imagFilename
auxL=get(handles.popupmenuPyramidLevel,'Value');
if(auxL>1)
    for ii=1:2*maxSize+1
        qq=ii+centerImgIdx-maxSize-1;
        if(qq>0 && qq<=numFrames)
            imgFilename{qq}=[imgFilename{qq}(1:end-4) '_GaussianPyramidLevel' num2str(auxL-1) fileExtension];
        end
    end
end


%{
%to use multiple cores (dangerous with memory)
if(matlabpool('size')~=8)
    if(matlabpool('size')>0)
        matlabpool close;
    end
    matlabpool(8);
end

stackAux=cell(2*handles.stackMaxSizeCache+1,1);
parfor ii=1:2*maxSize+1
    %stackAux{ii-centerImgIdx+maxSize+1}=handles.stackFilename{ii};
    qq=ii+centerImgIdx-maxSize-1;           
    if(qq>0 && qq<=numFrames)
        stackAux{ii}=load_v3d_raw_img_file(imgFilename{qq});
        %stackGlobal{ii}=load_v3d_raw_img_file(imgFilename{qq});
    %else
        %stackGlobal{ii}=[];
    end
end
stackGlobal=stackAux;
clear stackAux;
%}

%for single processor



for ii=1:2*maxSize+1
    %stackAux{ii-centerImgIdx+maxSize+1}=handles.stackFilename{ii};
    qq=ii+centerImgIdx-maxSize-1;           
    if(qq>0 && qq<=numFrames)
        switch( readType )
            case 0 
                stackGlobal{ii}=load_v3d_raw_img_file(imgFilename{qq});
                
                %disp('======Flipping y dimension to match TIFF matlab vs Gene lib========')%uncomment this if original image were tiff (it is not flipping x,y but doing Ny-y
                %stackGlobal{ii} = flipdim( stackGlobal{ii}, 2);
            case 1
                %disp('======PERMUTING JP2 FILES!!!===========')
                %stackGlobal{ii} = permute( readJPEG2000stack(imgFilename{qq},8),[2 1 3]);%to load the same way as load_v3d_raw_img_file
                disp('=========NOT PERMUTING JP2 FILES!!!===============')
                stackGlobal{ii} = readJPEG2000stack(imgFilename{qq},8);%if C code works directly on jp2 images flip is not necessary
            case 2
                stackGlobal{ii} = readKLBstack(imgFilename{qq});
            otherwise
                error 'Data type not readable'
        end
    else
        stackGlobal{ii}=[];
    end
end
%-------------------------------------------
handles.stackIniCache=centerImgIdx-handles.stackMaxSizeCache;
handles.stackFinCache=centerImgIdx+handles.stackMaxSizeCache;

set(handles.messageText,'String','All images loaded succesfully');
guidata(hObject,handles);
drawnow();

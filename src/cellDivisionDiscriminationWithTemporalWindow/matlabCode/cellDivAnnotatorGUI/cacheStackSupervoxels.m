function svFilename = cacheStackSupervoxels(handles,centerImgIdx, pathname)

global stackGlobal;
global svStructGlobal;


maxSize=handles.stackMaxSizeCache;


D = dir([pathname '*.svb']);
svFilename = cell(length(D) + handles.frameIni,1);
for ii = 1:length(D)
    svFilename{ii + handles.frameIni} = [pathname D(ii).name];
end
numFrames=length(svFilename);

addpath([fileparts( mfilename('fullpath') ) filesep 'readTGMM_XMLoutput' filesep 'readTGMM_XMLoutput'])
for ii=1:2*maxSize+1
    qq = ii+centerImgIdx-maxSize-1;           
    if(qq>handles.frameIni && qq<=numFrames)
        [svStructGlobal{ii}, sizeIm] = readListSupervoxelsFromBinaryFile( svFilename{qq} );
        
        if( isempty(sizeIm) == false )
            if( norm(sizeIm - size ( stackGlobal{ (length(stackGlobal) +1 )/2 } ) ) ~= 0 )
                sizeIm
                stackGlobal
                error 'Dimensions do not match'
            end
        else
            warning 'Supervoxels binary files not found. You cannot display them'
        end
    else
        svStructGlobal{ii}=[];
    end
end
rmpath([fileparts( mfilename('fullpath') ) filesep 'readTGMM_XMLoutput' filesep 'readTGMM_XMLoutput'])
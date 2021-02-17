%trackingCell{ii} = N x 18 array containing all the information related to i-th frame
% 
% Columns     1:10:   typical CATMAID array style
% Columns     11:16:  precision matrix W
% Columns     17:18:  childrenIdx. -1 indicates no children
% Columns     19:30:  svIdx for each cell. -1 indicate no sv
% Columns     31:31:  probability of Celldivision



%[trackingCell, neighCell];
function cacheTGMMxmlSolution(pathname, frame, frameIni, numNeigh)


global blobStructGlobal;
global neighGlobal

%check if the frame has already been loaded
if( frame == length(blobStructGlobal) - 1)%last frame
    if( isempty( blobStructGlobal{frame+1} ) == false )
        return;
    end    
elseif( isempty( blobStructGlobal{frame+1} ) == false && isempty( blobStructGlobal{frame+2} ) == false )%we need children in our display
    return;
end

addpath([fileparts( mfilename('fullpath') ) filesep 'readTGMM_XMLoutput' filesep 'readTGMM_XMLoutput'])

frameIni = max([frame - 1, frameIni]);%in case we need to set parent relationship
frameEnd = min([frame + 2, length(blobStructGlobal)-1]);%in case we need to set children relationship
basename = [pathname 'GMEMfinalResult_frame'];



%read frame information
[trackingMatrix, svIdxCell, precisionMatrix, probCellDivision]= parseMixtureGaussiansXml2trackingMatrixCATMAIDformat(basename,frameIni,frameEnd);

rmpath([fileparts( mfilename('fullpath') ) filesep 'readTGMM_XMLoutput' filesep 'readTGMM_XMLoutput'])

%parse svIdxCell into a constat array
maxSvPerNucleus = 12;%at the most a nuclei can have 12 supervoxels

%from http://stackoverflow.com/questions/6210495/how-to-combine-vectors-of-different-length-in-a-cell-array-into-matrix-in-matlab
maxLength = max(cellfun(@(x)numel(x),svIdxCell));
svIdxMatrix = cell2mat(cellfun(@(x)cat(2,x,-ones(1,maxLength-length(x))),svIdxCell,'UniformOutput',false));
if( size(svIdxMatrix,2) > maxSvPerNucleus )
   svIdxMatrix = svIdxMatrix(:,1: maxSvPerNucleus);
elseif( size(svIdxMatrix,2) < maxSvPerNucleus )
    svIdxMatrix = [svIdxMatrix, -ones( size(svIdxMatrix,1), maxSvPerNucleus - size(svIdxMatrix,2) )];
end


%organize elements into one cell per time point
tIni = min(trackingMatrix(:,8));
tEnd = max(trackingMatrix(:,8));

[p, ~] = hist(trackingMatrix(:,8),[tIni:tEnd]);
trackingCell = mat2cell(single([trackingMatrix precisionMatrix -ones(size(precisionMatrix,1),2) svIdxMatrix probCellDivision]),p);%we do not need double precision
%update parentId relationship
ii = 1;
treenodeIdMap = containers.Map(trackingCell{ii}(:,1),[1:size(trackingCell{ii},1)]');

for ii = 2:length(trackingCell)
    parIdCell = num2cell( trackingCell{ii}(:,7) );
    
    mask = isKey(treenodeIdMap, parIdCell);
    
    trackingCell{ii}(~mask,7) = -1;%no parent
    trackingCell{ii}(mask,7) = cell2mat(values(treenodeIdMap,parIdCell(mask)))-1;%C-indexing to make it easier to port the code from blobStruct to CATMAID array style
    
    
    for jj = 1:size(trackingCell{ii},1)
        parIdx = trackingCell{ii}(jj,7) + 1;
        if( parIdx > 0 )
           if( trackingCell{ii-1}(parIdx,17) == -1 )%no children yet
               trackingCell{ii-1}(parIdx,17) = jj - 1;%C-indexing
           else%cell division
               trackingCell{ii-1}(parIdx,18) = jj - 1;%C-indexing
           end
        end
    end
    
    %update hash for next round
    treenodeIdMap = containers.Map(trackingCell{ii}(:,1),[1:size(trackingCell{ii},1)]');
end

%copy final structure
if( frame == frameIni )
    blobStructGlobal{frame + 1} = trackingCell{1};
    blobStructGlobal{frame + 2} = trackingCell{2};
else
    blobStructGlobal{frame + 1} = trackingCell{2};
    blobStructGlobal{frame + 2} = trackingCell{3};
end
%calculate nearest neighbors
if(numNeigh > 0 )
   eror 'TODO: nearest neighbor not implemented yet' 
end
    
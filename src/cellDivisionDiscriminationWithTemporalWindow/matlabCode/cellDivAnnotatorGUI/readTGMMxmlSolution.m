%trackingCell{ii} = N x 18 array containing all the information related to i-th frame
% 
% Columns     1:10:   typical CATMAID array style
% Columns     11:16:  precision matrix W
% Columns     17:18:  childrenIdx. -1 indicates no children
% Columns     19:30:  svIdx for each cell. -1 indicate no sv



function [trackingCell, neighCell] = readTGMMxmlSolution(pathname, numNeigh)

addpath([fileparts( mfilename('fullpath') ) filesep 'readTGMM_XMLoutput' filesep 'readTGMM_XMLoutput'])

D = dir([pathname 'GMEMfinalResult_frame*.xml']);

frameIni = str2double(D(1).name(end-7:end-4));
frameEnd = str2double(D(end).name(end-7:end-4));
basename = [pathname 'GMEMfinalResult_frame'];

[trackingMatrix, svIdxCell, precisionMatrix]= parseMixtureGaussiansXml2trackingMatrixCATMAIDformat(basename,frameIni,frameEnd);

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
trackingCell = mat2cell(single([trackingMatrix precisionMatrix -ones(size(precisionMatrix,1),2) svIdxMatrix]),p);%we do not need double precision
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

%calculate nearest neighbors
neighCell = cell(size(trackingCell,1));
if(numNeigh > 0 )
   eror 'TODO: nearest neighbor not implemented yet' 
end
    
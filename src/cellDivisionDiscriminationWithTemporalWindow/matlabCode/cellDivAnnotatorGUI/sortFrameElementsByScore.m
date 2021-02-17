%sorts elements of a frame by size

function [sortIdxList, aux]=sortFrameElementsByScore(blobStructFrame)



vol=zeros(length(blobStructFrame),1);
for ii=1:length(blobStructFrame)
    
    if(isempty(blobStructFrame(ii).frame))
        vol(ii)=-1e32;
    else
        vol(ii) = blobStructFrame(ii).solutions(1).score;        
    end
end

[aux sortIdxList]=sort(vol,'descend');%returns 

pos = find( aux == aux(1) );%randomize top score
if( length(pos) > 1)
    pp = randperm( length(pos) );
    
    aux( pos ) = aux( pos ( pp ) );
    sortIdxList( pos ) = sortIdxList( pos ( pp ) );
end

%{
%%for structure learning annotions to find very similar ones

vol= 1e32 * ones(length(blobStructFrame),1);
for ii=1:10:length(blobStructFrame)-1
    
    if(isempty(blobStructFrame(ii).frame) || isempty(blobStructFrame(ii+1).frame) || blobStructFrame(ii).solutions(1).score > 1e30 || blobStructFrame(ii).solutions(1).score < 0 )
        vol(ii)= 1e32;
    else
        vol(ii) = abs(blobStructFrame(ii+1).solutions(1).score - blobStructFrame(ii).solutions(1).score);        
    end
end
vol(end) = 1e32;
[aux sortIdxList]=sort(vol,'ascend');%returns 

%}
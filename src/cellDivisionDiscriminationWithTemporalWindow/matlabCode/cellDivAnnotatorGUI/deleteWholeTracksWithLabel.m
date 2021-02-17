function blobStruct=deleteWholeTracksWithLabel(blobStruct,labelList)

dims=length(blobStruct(1,1).center);

count=0;
for ii=1:size(blobStruct,1)
    for jj=1:size(blobStruct,2)
        if(isempty(blobStruct(ii,jj).frame))
            break;
        end
        label=blobStruct(ii,jj).solutions(1).label;
        if(sum(label==labelList)>0)%kill cell
            blobStruct(ii,jj).solutions(1).label=0;
            blobStruct(ii,jj).solutions(1).parentIdx=[uint32(2^32-1) 0];
            blobStruct(ii,jj).solutions(1).childrenIdx=[];
            blobStruct(ii,jj).center=repmat(-1e32,[1 dims]);
            blobStruct(ii,jj).surface.coeffs(7:9)=repmat(-1e32,[1 dims]);
            count=count+1;
        end
    end
end

disp(['Deleted ' num2str(count) ' cells'])



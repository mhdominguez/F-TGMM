function numObj = parseMixtureGaussiansXml2MCMCBinary_MultiFrame(basename,imgFilenamePattern,frameIni,frameEnd, maxKnn, fixedTimePoint)

if(exist('maxKnn','var') == 0)
    maxKnn=10;
end


if(exist('fixedTimePoint','var') == 0)
    fixedTimePoint = -1; %to use a unique time point for all TGMM files. Especial case to debug hierarchical segmentation
end

fout=fopen([basename '_stackMCMC_TM' num2str(frameIni,'%.4d') '_' num2str(frameEnd,'%.4d')  '.bin'],'wb');


maxLabel=0;
numObj = zeros(length(frameIni:frameEnd),1);

%write indetifier char
fwrite(fout,'blobstruct','char');
%write number of frames
fwrite(fout,int32(frameEnd-frameIni+1),'int32');
%write number of neighbors
fwrite(fout,int32(maxKnn),'int32');


for frame=frameIni:frameEnd
    
    if( fixedTimePoint < 0 )
        frameIm = frame;
    else
        frameIm = fixedTimePoint;
    end
    
    
    imgFilename = recoverFilenameFromPattern(imgFilenamePattern, frame);
    if(exist(imgFilename,'file')==0)
        error(['Image filename ' imgFilename ' does not exist']);
    end
    
    if(frame==frameIni)
        obj=readXMLmixtureGaussians([basename num2str(frame,'%.4d') '.xml']);
    
        %write dimensions only once
        fwrite(fout, obj(1).dims,'int32');
        %write scale only once
        fwrite(fout, single(obj(1).scale), 'float32');%scale
        
        %-----------clean dead cells----------------
        %mapId(obj(kk).id+1)=obj(kk).newId; if obj(kk) is
        %dead->mapId(obj(kk).id)=-1
        count=0;
        erase=[];
        mapId=-ones(length(obj),1);
        for ii=1:length(obj)
            if(obj(ii).m<-1e31)
                erase=[erase;ii];
            else
                mapId(ii)=count;%C indexing
                count=count+1;
            end
        end        
        if(~isempty(erase))
            obj(erase)=[];
            %mapId(erase)=[];
            display(['Deleted ' num2str(length(erase)) ' cells in frame ' num2str(frame)])
        end
        %----------------------------------------------
        mapIdPar=[];
    else
        obj=objCh;
        mapIdPar=mapId;
        mapId=mapIdCh;
    end
    
    if(frame~=frameEnd)
        objCh=readXMLmixtureGaussians([basename num2str(frame+1,'%.4d') '.xml']);
    else
        objCh=[];
    end
    
    %------------clean dead cells--------------------------
    count=0;
    erase=[];
    mapIdCh=-ones(length(objCh),1);
    for ii=1:length(objCh)
        if(objCh(ii).m<-1e31)
            erase=[erase;ii];
        else
            mapIdCh(ii)=count;%C indexing
            count=count+1;
        end
    end
    if(~isempty(erase))
        objCh(erase)=[];
        %mapIdCh(erase)=[];
        display(['Deleted ' num2str(length(erase)) ' cells in frame ' num2str(frame+1)])
    end
    %----------------------------------------------------
    
    
    %figure out parents
    objCh_parents=zeros(length(objCh),1);
    for kk=1:length(objCh_parents)
       objCh_parents(kk)=objCh(kk).parent; 
    end
    
    %write imgfilenae   
    fwrite(fout,[imgFilename '*'],'char');% chacter * is the ending character
    
    %write supervoxel filename
    svFile = [basename num2str(frame,'%.4d') '.svb'];
    fwrite(fout,[svFile '*'],'char');% chacter * is the ending character
    
    %perform nearest neighbors to facilitate display
    mm=zeros(length(obj),3);
    for kk=1:size(mm,1)
        mm(kk,:)=obj(kk).m.*obj(1).scale;
    end
    idx=knnsearch(mm,mm,'k',maxKnn+1);%maxKnn+1 neighbors per Gaussian
    idx(:,1)=[];%first column is itself
    
    if( size(idx, 2) < maxKnn )%less points than nearest neighbors in this frame
        idx = [ idx, ones(size(idx,1), maxKnn - size(idx,2) ) ];
    end
    
    %write number of objects to read
    fwrite(fout, length(obj),'int32');
    %write frame
    fwrite(fout, frame-frameIni, 'int32');
    
    %process each object
    for ii=1:length(obj)
        blob=obj(ii);
        maxLabel=max(maxLabel,blob.lineage+1);
        if(isempty(blob.splitScore)) 
            blob.splitScore=-1e32;
        end;
        %blob=getBlobStructure(aux(:,1),size(img,2)-1-aux(:,2),aux(:,3),img,scale);%since index goes from [0,N-1] in order to flip we have to do N-1-y
        
        %write XML blob
        %{
        fprintf(fout,'<Blob id="%d" dims="%d" center="%f %f %f" scale="%.2f %.2f %.2f" frame="%d" intensity="%.2f" neigh="',mapId(blob.id+1),blob.dims,blob.m(1),blob.m(2),blob.m(3),...
            blob.scale(1),blob.scale(2),blob.scale(3),frame-frameIni,max(blob.alpha-blob.alphaPrior,0));
        for kk=1:maxKnn
            fprintf(fout,'%d %d ',frame-frameIni,mapId(obj(idx(ii,kk)).id+1));%C-indexing
        end
        fprintf(fout,'">\n');
        %}
        fwrite(fout, mapId(blob.id+1), 'uint32');%blobid
        fwrite(fout, single(blob.m), 'float32');%center
        fwrite(fout, single(max(blob.alpha-blob.alphaPrior,0)), 'float32');%intensity
        fwrite(fout, single(blob.betaPrior), 'float32');%used as a wildcard for different measuemrents: betaPrior has probBAckground
        for kk=1:maxKnn
            fwrite(fout,[frame-frameIni,mapId(obj(idx(ii,kk)).id+1)], 'uint32');%C-indexing
        end
        
        %{
        %write supervoxel idx
        fprintf(fout,'" svIdx=" ');
        fprintf(fout,'%d ', blob.svIdx);%C-indexing
        %}
        fwrite(fout, int32(length(blob.svIdx)), 'int32');
        if( isempty( blob.svIdx ) == false )
            fwrite(fout, int32(blob.svIdx), 'int32');
        end
        
        %{
        fprintf(fout,'<Surface name="Ellipsoid" id="1" numCoeffs="%d"\n',blob.dims+(blob.dims*(blob.dims+1))/2);
        fprintf(fout,'coeffs="');
        fprintf(fout,'%f ',[blob.nu*[blob.W(1,1) blob.W(1,2) blob.W(1,3) blob.W(2,2) blob.W(2,3) blob.W(3,3)] blob.m]);
        fprintf(fout,'" covarianceMatrixSize="%d" > </Surface>\n',blob.dims);
        %}
        %write num coeffs
        fwrite(fout, blob.dims+(blob.dims*(blob.dims+1))/2, 'int32');
        %write coeffs 
        fwrite(fout, single([blob.nu*[blob.W(1,1) blob.W(1,2) blob.W(1,3) blob.W(2,2) blob.W(2,3) blob.W(3,3)] blob.m]), 'float32');
        
        ch=find(objCh_parents==blob.id);
        if(isempty(ch) || frame==frameEnd)%no children            
            if(frame==frameIni || blob.parent < 0 )%no parents
                %fprintf(fout,'<BlobSolution score="%g" label="%d" parentIdx="%u %u"> </BlobSolution>\n',blob.splitScore,blob.lineage+1,uint32(2^32),uint32(0));%current solution
                %fprintf(fout,'<BlobSolution score="%g" label="%d" parentIdx="%u %u"> </BlobSolution>\n',blob.splitScore,blob.lineage+1,uint32(2^32),uint32(0));%minimum energy solution
                
                %WE WRITE ONLY ONCE!!
                fwrite(fout, single(blob.splitScore), 'float32');
                fwrite(fout, blob.lineage+1, 'int32');
                fwrite(fout, [uint32(2^32),uint32(0)],'uint32');
                fwrite(fout,0,'int32');%number of children
            else%parents
                %fprintf(fout,'<BlobSolution score="%g" label="%d" parentIdx="%u %u"> </BlobSolution>\n',blob.splitScore,blob.lineage+1,uint32(frame-frameIni-1),uint32(mapIdPar(blob.parent+1)));%current solution
                %fprintf(fout,'<BlobSolution score="%g" label="%d" parentIdx="%u %u"> </BlobSolution>\n',blob.splitScore,blob.lineage+1,uint32(frame-frameIni-1),uint32(mapIdPar(blob.parent+1)));%minimum energy solution
                %WE WRITE ONLY ONCE!!
                fwrite(fout, single(blob.splitScore), 'float32');
                fwrite(fout, blob.lineage+1, 'int32');
                fwrite(fout, [uint32(frame-frameIni-1),uint32(mapIdPar(blob.parent+1))],'uint32');
                fwrite(fout,0,'int32');%number of children
                
            end
        else
            if(frame==frameIni || blob.parent < 0 )%no parents
               
                %WE WRITE ONLY ONCE!!
                fwrite(fout, single(blob.splitScore), 'float32');
                fwrite(fout, blob.lineage+1, 'int32');
                fwrite(fout, [uint32(2^32),uint32(0)],'uint32');
                fwrite(fout,length(ch),'int32');%number of children
                for cc=1:length(ch)
                    fwrite(fout, [uint32(frame-frameIni+1),uint32(mapIdCh(objCh(ch(cc)).id+1))], 'uint32');
                end
                
            else%parents
                
                
                %WE WRITE ONLY ONCE!!
                fwrite(fout, single(blob.splitScore), 'float32');
                fwrite(fout, blob.lineage+1, 'int32');
                fwrite(fout, [uint32(frame-frameIni-1),uint32(mapIdPar(blob.parent+1))],'uint32');
                fwrite(fout,length(ch),'int32');%number of children
                for cc=1:length(ch)
                    fwrite(fout, [uint32(frame-frameIni+1),uint32(mapIdCh(objCh(ch(cc)).id+1))], 'uint32');
                end
                
            end
        end
        %fprintf(fout,'</Blob>\n');
        
        if(blob.m(1)>=0)%object not dead
            numObj(frame-frameIni+1)=numObj(frame-frameIni+1)+1;
        end
    end
    
    %write footer for frame
   % fprintf(fout,'</Frame>');
   
    
end

%write footer
%fprintf(fout,'</Stack>\n');
%fprintf(fout,'</document>\n');

fclose(fout);

disp(['PLEASE open the file ' basename '_stackMCMC.bin' ' and change maxLabels attribute to ' num2str(maxLabel+1)]);

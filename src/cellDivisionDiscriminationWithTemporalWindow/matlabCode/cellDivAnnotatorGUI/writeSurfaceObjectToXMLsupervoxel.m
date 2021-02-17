%we do not include xml footer and header, so you have to add it yourself
%before hand

function writeSurfaceObjectToXMLsupervoxel(obj,fout)

if(isempty(obj)) return;end;

numCoeffs=length(obj(1).coeffs);

for kk=1:length(obj)
    blob=obj(kk);
    fprintf(fout,'<Surface name="Ellipsoid" id="1" numCoeffs="%d" ',numCoeffs);
    fprintf(fout,'coeffs="');
    fprintf(fout,'%f ',blob.coeffs);
    
    if (isfield(blob,'intensity') == true)
        fprintf(fout,'" intensity="%g',blob.intensity);
    end
    
    fprintf(fout,'" covarianceMatrixSize="%d" \n',blob.covarianceMatrixSize);
    fprintf(fout,'svFilename="%s" svIdx="',blob.svFilename);
    
    fprintf(fout,'%d ',blob.svIdx');
    
    fprintf(fout,'"\n imFilename="%s" class="%s"></Surface>\n',blob.imFilename,blob.class);    
end
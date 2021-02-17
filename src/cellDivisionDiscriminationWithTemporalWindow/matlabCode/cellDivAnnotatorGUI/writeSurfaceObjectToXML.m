%we do not include xml footer and header, so you have to add it yourself
%before hand

function writeSurfaceObjectToXML(obj,fout)

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
    fprintf(fout,'imFilename="%s" class="%s"></Surface>\n',blob.imFilename,blob.class);    
end
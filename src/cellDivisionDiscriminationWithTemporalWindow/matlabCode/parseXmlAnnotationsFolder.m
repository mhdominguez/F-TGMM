%we assume surface is 3D
function [objFinal, ff] = parseXmlAnnotationsFolder(folderPath)



%get list of xml files
files=dir([folderPath filesep '*.xml']);


N=length(files);
objFinal = [];
for kk=1:N
    %read xml file
    xDoc = xmlread([folderPath filesep files(kk).name]);
    
    
    NodeList = xDoc.getElementsByTagName('Surface');
    numG=NodeList.getLength();
    
    obj(numG,1).imFilename=' ';%preallocate
    
    
    for ii=1:numG
        fstNode = NodeList.item(ii-1);%contains the ii-th Surface
        attrs = fstNode.getAttributes();
        
        if(getNumericalAttribute(attrs,'covarianceMatrixSize')~=3)
            error 'Dimensions is not 3D';
        end
        aux=getNumericalAttribute(attrs,'coeffs');
        obj(ii).m=aux(7:9);%center of the ellipsoid 
        obj(ii).W=[aux(1:3); aux(2) aux(4:5); aux(3) aux(5) aux(6)];
        obj(ii).imFilename=char(attrs.getNamedItem('imFilename').getValue);
        obj(ii).class=char(attrs.getNamedItem('class').getValue);        
    end
    objFinal = vertcat(objFinal,obj);
    clear obj;
end


%calculate histogram for different classes
ff = zeros(length(objFinal),1);
for kk=1:length(ff)
   switch objFinal(kk).class
       
       case 'cellDivisionWrong'
           ff(kk)=0;
       case 'cellDivisionCorrect'
           ff(kk)=1;
   end
end



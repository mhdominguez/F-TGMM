%Fernando Amat July 16th 2010
%Reads a multiple-image TIFF stac
function [stack info]=readTIFFstack(filename)

info=imfinfo(filename);

N=length(info);

im=imread(filename,'Info',info,'Index',1);

if(info(1).SamplesPerPixel==1)
[W H]=size(im);
else
    [W H D]=size(im);
end

type_=class(im);

if(info(1).SamplesPerPixel==1)
    stack=zeros(W,H,N,type_);
    stack(:,:,1)=im;
    for kk=2:N
        stack(:,:,kk)=imread(filename,'Info',info,'Index',kk);
    end
else
    stack=zeros(W,H,N,info(1).SamplesPerPixel,type_);%RGB
    for ii=1:info(1).SamplesPerPixel
        stack(:,:,1,ii)=im(:,:,ii);
    end
    for kk=2:N
        for ii=1:info(kk).SamplesPerPixel
            aux=imread(filename,'Info',info,'Index',kk);
            stack(:,:,kk,ii)=aux(:,:,ii);
        end
    end
end



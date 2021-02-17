%substitute to avoid installed V3D toolbox and preserve old code
function im=load_v3d_raw_img_file(imgFilename)

[im, info]=readTIFFstack(imgFilename);


if(info(1).SamplesPerPixel==1)%gray image
   
    if(info(1).BitDepth==8)%UINT8: y is flipped
        if(length(size(im))==2)%2D image with RGB
            im=permute(im,[2 1]);
        else
            im=permute(im,[2 1 3:length(size(im))]);
        end
        im=flipdim(im,2);
    elseif(info(1).BitDepth==16)%UINT16
        if(length(size(im))==2)%2D image with RGB
            im=permute(im,[2 1]);
        else
            im=permute(im,[2 1 3:length(size(im))]);
        end
    else
        %at this point I do not knwo what to
    end
    
else%RGB image assuming 8 bits per channel

    im=permute(im,[2 1 3:length(size(im))]);
    
    
end


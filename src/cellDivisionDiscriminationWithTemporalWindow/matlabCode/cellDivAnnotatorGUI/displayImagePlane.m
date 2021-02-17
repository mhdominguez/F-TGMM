%displays image plane on top of cell lineages

%INPUT: imgzposition: desired z position of the image plane.

function hsurf=displayImagePlane(imgzposition,im,hAxis,hSurfIm,rectCoord)

% scale image between [0, 255] in order to use a custom color map for it.
if(isa(im,'uint8')==0)
    im=single(im);
    im=im-min(im(:));
    im=uint8(255*im/max(im(:)));
end
% convert the image to a true color image with the jet colormap.
colorimg = ind2rgb(im',gray(256));

% plot the image plane using surf.
if(~isempty(hSurfIm))
    delete(hSurfIm);
end

hold(hAxis,'on');
hsurf=surf(hAxis,[rectCoord(1) rectCoord(2)],[rectCoord(3) rectCoord(4)]-1,repmat(imgzposition-1, [2 2]),colorimg,'facecolor','texture');%to adjust rectCoord and imgZposition to C-indexing
hold(hAxis,'off');
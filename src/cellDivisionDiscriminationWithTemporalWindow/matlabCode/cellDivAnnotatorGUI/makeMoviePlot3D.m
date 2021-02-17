%makes a movie roating around any angle defined by view with the current
%figure

%example from http://www.mathworks.com/help/techdoc/ref/movie.html
function F=makeMoviePlot3D(angleStep,numFrames)

%axis tight
set(gca,'nextplot','replacechildren');
% Record the movie
[az el]=view;
for j = 1:numFrames 
    %surf(sin(2*pi*j/20)*Z,Z)
    
    %F(j) = getframe(gcf,get(gcf,'position'));
    F(j) = getframe(gcf);
    view([az+angleStep*j el]);
end


%{
%to visualize movie and save it
[h, w, p] = size(F(1).cdata);  % use 1st frame to get dimensions
hf = figure;
% resize figure based on frame's w x h, and place at (150, 150)
set(hf, 'position', [150 150 w h]);
axis off
% tell movie command to place frames at bottom left
movie(hf,F,4,5,[0 0 0 0]);
%}

%to save the movie
%movie2avi(F,'/Users/amatf/TrackingNuclei/reports/videos/solutionComparison.avi','FPS',10,'Compression','None');
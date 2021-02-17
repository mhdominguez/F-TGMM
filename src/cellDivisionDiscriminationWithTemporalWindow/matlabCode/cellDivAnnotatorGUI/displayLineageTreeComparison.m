%displays different lineages from different solutions
function displayLineageTreeComparison(blobStructCell,frame,blob,numSolution,anisotropy)

N=length(blobStructCell);
M=length(frame);

h=figure;
countC=0;
for ii=1:M
    for kk=1:N
        displayLineageTree(blobStructCell{kk},frame(ii),blob(ii),numSolution(kk),anisotropy,gca,countC,(kk-1)*[4 4 0]);
        countC=countC+1;
    end
end
colormap(lines(N*M));
colorbar
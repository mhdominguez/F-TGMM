par = [30 2512] + 1
ch = blobStructGlobal(par(1) , par(2)).solutions.childrenIdx;

%go to root of the branch
while( par(1) < 4e9 && length(ch) <= 2 )     
   parOld = par;
   par = blobStructGlobal(par(1) , par(2)).solutions.parentIdx + 1;
   if( par(1) > 4e9 ) break;end;
   ch = blobStructGlobal(par(1) , par(2)).solutions.childrenIdx + 1;
end

if( length(ch) > 2 || par(1) > 4e9)%move forward once
    par = parOld;
end

%%
%traverse the branch forward
rootBranch = par;
ch = rootBranch;
bb = [];
while( length(ch) == 2 )
    bb = [ bb; [ch(1)-1 ch(2)-1 blobStructGlobal(ch(1) , ch(2)).radius ]];
    ch = blobStructGlobal(ch(1) , ch(2)).solutions.childrenIdx + 1;
end
function [sortIdxList, aux]=sortFrameElementsByProbCellDivision(blobStructFrame,blobStructFrameCh)

numObj=size(blobStructFrame,1);

vol=-ones(numObj,1);
for ii=1:numObj
    
    
    
    
    
    if(blobStructFrame(ii,3)<-1e30)%dead cell
        vol(ii)=-1e32;
    elseif(blobStructFrame(ii,18) >= 0)%cell division
        %make sure both children are alive                
        ch=blobStructFrame(ii,[17 18]);
        %two children
        ch = [-1, ch(1), -1, ch(2)];        
        
        if(blobStructFrameCh(ch(2)+1, 3)>-10 && blobStructFrameCh(ch(4)+1, 3)>-10)
            
            vol(ii) = blobStructFrame(ii,31);%probability of cell division
        else
            vol(ii)=-1e32;
        end
    end
    
end

[aux, sortIdxList]=sort(vol,'descend');%returns

% % % % %randomized the order of teh cell division
% % % % pos=find(aux==1);
% % % % sortIdxList(pos)=sortIdxList(pos(randperm(length(pos))));
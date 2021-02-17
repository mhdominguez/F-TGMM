function [sortIdxList, aux]=sortFrameElementsByCellDivision(blobStructFrame,blobStructFrameCh)

numObj=size(blobStructFrame,1);

vol=-ones(numObj,1);
for ii=1:numObj
    
    
    
    
    
    if(blobStructFrame(ii,3)<-1e30)%dead cell
        vol(ii)=-1e32;
    elseif(blobStructFrame(ii,18) >= 0)
        %make sure both children are alive                
        ch=blobStructFrame(ii,[17 18]);
        %two children
        ch = [-1, ch(1), -1, ch(2)];        
        
        if(blobStructFrameCh(ch(2)+1, 3)>-10 && blobStructFrameCh(ch(4)+1, 3)>-10)
            %{
                %vol(ii)=1;
                %calculate Mahalanobis distance between daughhters and
                %mother cell
                ccMom = blobStructFrame(ii).center;
                
                WD = blobStructFrameCh(ch(2)+1).surface.coeffs(1:6);
                WD1 = [WD(1:3); WD(2) WD(4) WD(5);WD(3) WD(5) WD(6)];
                ccD1 = blobStructFrameCh(ch(2)+1).surface.coeffs(7:9);
                
                WD = blobStructFrameCh(ch(4)+1).surface.coeffs(1:6);
                WD2 = [WD(1:3); WD(2) WD(4) WD(5);WD(3) WD(5) WD(6)];
                ccD2 = blobStructFrameCh(ch(4)+1).surface.coeffs(7:9);
                
                p1 = (ccD1-ccMom) * WD1 * (ccD1-ccMom)';
                p2 = (ccD2-ccMom) * WD2 * (ccD2-ccMom)';
                
                logitS = -0.5* ( -log( det(WD1) ) + p1 + log( det(WD2) ) -p2 );
                
                %apply absolute value since we are looking for cases were
                %it is not a clear cut
                logitS = abs(logitS);
                
                vol(ii) = 1.0/(logitS + eps);%higher ones should be the real (or more symmetric divisions)
                %TODO: add mass eveness here and weight the score somehow
            %}
            
            
            %check how close they are to mid point between daughters
            ccMom = blobStructFrame(ii,3:5);
            ccD1 = blobStructFrameCh(ch(2)+1,3:5);
            ccD2 = blobStructFrameCh(ch(4)+1,3:5);
            pp = 0.5 * (ccD1 + ccD2);
            v1 = ccD1 - ccD2;%to construct covariance (we penalize less moving along the axis than perpendicular to it
            v1 = v1 / norm(v1);
            if( abs(v1(1)) < 0.98 )
                v2 = [1 0 0 ];%make v1 and v2 ar eindependent
            else
                v2 = [0 1 0];
            end
            v2 = v2 - dot(v1,v2) * v1;
            v2 = v2 / norm(v2);
            v3 = cross(v1,v2);
            W = [v1' v2' v3'] * diag(1./[2, 1 , 1]) * [v1;v2;v3];
            dd = (ccMom - pp) * W * (ccMom - pp)';
            vol(ii) = 1.0 / (dd + eps);
        else
            vol(ii)=-1e32;
        end
    end
    
end

[aux, sortIdxList]=sort(vol,'descend');%returns

% % % % %randomized the order of teh cell division
% % % % pos=find(aux==1);
% % % % sortIdxList(pos)=sortIdxList(pos(randperm(length(pos))));
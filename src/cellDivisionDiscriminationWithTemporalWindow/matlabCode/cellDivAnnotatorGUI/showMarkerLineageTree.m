%========================================================================
%display marker to show where we are in the lineage tree
function hDot=showMarkerLineageTree(blob,h,hDot,zLevel)

disp('WARNING:hDot=showMarkerLineageTree not working properly. Deactivated until it is fixed')
%{
zStep=30;%this value has to match displayLineageTree.m
offset=0;%this value has to match displayLineageTree.m

if(isempty(zLevel))
    zLevel=0;
end

if(~isempty(hDot))
    delete(hDot);
end

sigma1=zeros(3);
sigma1(1,1)=blob.surface.coeffs(1);sigma1(2,1)=blob.surface.coeffs(2);sigma1(3,1)=blob.surface.coeffs(3);
sigma1(2,2)=blob.surface.coeffs(4);sigma1(3,2)=blob.surface.coeffs(5);sigma1(3,3)=blob.surface.coeffs(6);
sigma1(1,2)=sigma1(2,1);sigma1(1,3)=sigma1(3,1);sigma1(2,3)=sigma1(3,2);

[V D]=eig(sigma1);

e1=V(:,1)/sqrt(D(1));

hold(h,'on');hDot=plot3(h,blob.center(1)+e1(1),blob.center(2)+e1(2),blob.center(3)+offset-zLevel*zStep+e1(3),'r.','MarkerSize',100);hold(h,'off');
%WARNING: you need to update teh Zlevel outside depending if we increment
%or decrement in the lineage tree
%}
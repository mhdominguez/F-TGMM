function m = median3elem(a)
m = a(1+1);
if (a(0+1) <= a(1+1))%median of 3 elements: from http:%stackoverflow.com/questions/4793251/sorting-int-array-with-only-3-elements
    if (a(1+1) > a(2+1))
        
        if (a(0+1) < a(2+1))
            
            %temp = a(0+1);
            %a(0+1) = a(2+1);
            %a(2+1) = a(1+1);
            %a(1+1) = temp;
            m = a(2+1);
        else
            %temp = a(1+1);
            %a(1+1) = a(2+1);
            %a(2+1) = temp;
            m = a(0+1);
        end
    end
else
    if (a(0+1) > a(2+1))
        
        if (a(1+1) < a(2+1))
            
            %temp = a(0+1);
            %a(0+1) = a(1+1);
            %a(1+1) = a(2+1);
            %a(2+1) = temp;
            m = a(2+1);
        else
            %temp = a(0+1);
            %a(0+1) = a(2+1);
            %a(2+1) = temp;
            %auxNode->data = a(1+1); NO NEED TO REASSIGN
        end
    else
        %temp = a(0+1);
        %a(0+1) = a(1+1);
        %a(1+1) = temp;
        m = a(0+1);
    end
end
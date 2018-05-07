%Written by Neil McHenry
% April 16th 2018
%HW3p1.m

function [VS,R] = catalogsearch(range)

%Read Star Catalog
    load mappar;
%Sort the stars in Magnitude (ascending from brightest to dimmest)
    TOTAL = [MAGN;VS];
    [TOTAL_temp, INDEX] = sort(TOTAL(1,:));
    TOTALSORT = TOTAL(:,INDEX);
    Stars = TOTALSORT(1,:);
%% K-Vector to select stars between magnitude ranges
    %Define Variables
        %Number of Elements
        n = 10000;
        %Epsilon (scales by number of elements)
        eps =1/n;
        m = min(Stars)-eps;
        M = max(Stars)+eps;
        %Step Calculation
        p = ((M-m)/(n-1));
        %Offset calculation
        q = m-p;
        %K-Vector Calculation
        k=[];
        for i=1:n
            %Find # of elements below each step division
            numStars = numel(find(Stars<=(q+(p*(i)))));
            k(i)= numStars;
        end

%% n = 10,000 k-vector range searching
    %Specify random magnitude range
%         range = [-1.5,6];
        a = range(1);
        b = range(2);
    %Calculate the j values for the bottom and top of the range
        j_a = floor((a-q)/p);
        j_b = ceil((b-q)/p);
    %Calculate k start (k_a) and end (k_b)
        k_a = k(j_a)+1; while(Stars(k_a)<a && k_a<=k_b) k_a=k_a+1; end
        k_b = k(j_b); while(Stars(k_b)>b && k_b>=k_a) k_b=k_b-1; end
    %Save the range of values to the R Matrix
        R = TOTALSORT(1:4,k_a:k_b);
    
end
    

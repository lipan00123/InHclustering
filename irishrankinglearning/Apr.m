function [c1,c2]=Apr(mutualinfoest,n1,n2)
n=n1+n2;
val = zeros(n*(n-1)/2,1);
for i = 2:n,
    for j = 1:i-1,
        row = (i-1)*(i-2)/2+j;
        [tempsort, index] = sort(mutualinfoest(row,:));
        temprowlist = generate_rowlist(index(1:n1));
        for k = 1:length(temprowlist),
            for p = 1:n2,
                val(row) = val(row) + mutualinfoest(temprowlist(k),index(p+n1));
            end
        end
        temprowlist = generate_rowlist(index(n1+1:n1+n2));
        for k = 1:length(temprowlist),
            for p = 1:n1,
                val(row) = val(row) + mutualinfoest(temprowlist(k),index(p));
            end
        end
    end
end

[minval, minindex] = min(val);
[tempsort, index] = sort(mutualinfoest(minindex,:));
if length(find(tempsort==0))<=n1,
    c1 = sort(index(1:n1));
    c2 = sort(index(n1+1:n1+n2));
else
    temp = find(mutualinfoest(minindex,:)==0);
    
    c1 = temp(randperm(length(temp),n1));
    temp = ones(n,1);
    temp(c1) = 0 ;
    c2 = find(temp==1);
end
end
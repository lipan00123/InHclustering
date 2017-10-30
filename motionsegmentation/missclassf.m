
function [err,assignment]=missclassf(estlabel,label)
    c = max(label);
    cost = zeros(c,c);
    for i = 1:c,
        for j = 1:c,
            cost(i,j)=sum(estlabel(find(label==i))~=j);
        end
    end
    [assignment,err] = munkres(cost);
end
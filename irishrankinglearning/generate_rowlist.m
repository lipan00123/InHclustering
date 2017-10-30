function y=generate_rowlist(index)
    r = length(index);
    index=sort(index,'ascend');
    y = zeros(r*(r-1)/2,1);
    for i = 2:length(index),
        for j=1:i-1,
            y((i-1)*(i-2)/2+j)=(index(i)-1)*(index(i)-2)/2+index(j);
        end
    end
end
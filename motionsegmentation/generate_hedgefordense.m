function edgelist=generate_hedgefordense(label,p,edgenum)
%%%%%%%%% pick edge for dense %%%%%%%%%even in each data point
m = length(label);
K = max(label);
edgelist = zeros(edgenum,p);
a = unidrnd(m,[edgenum,1]);
for i = 1:edgenum,
    tempseg = find(label == label(a(i)));
    if length(tempseg)<p,
        edgelist(i,:) = datasample(1:m,p,'Replace',false);
    else
        edgelist(i,:) = datasample(tempseg,p,'Replace',false);
    end
end
end

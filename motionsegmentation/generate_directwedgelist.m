function directwedgelist=generate_directwedgelist(data, edgelist)


[edgenum,hedgesize]=size(edgelist);
[dim, m] = size(data);
directwedgelist = zeros(hedgesize+1,m-hedgesize,edgenum);
for i = 1:edgenum,
    U = zeros(dim,hedgesize-2,hedgesize);
    b = zeros(dim,hedgesize);
    A = data(:,edgelist(i,1:hedgesize));
    Apr = bsxfun(@minus, A(:,1:hedgesize-1) , A(:,hedgesize));
    dis = zeros(hedgesize+1,m-hedgesize);
    for j = 1:hedgesize-1,
        [U(:,:,j),D,V]=svd(Apr(:,[1:j-1 j+1:hedgesize-1]),'econ');
        b(:,j) = Apr(:,j)-U(:,:,j)*(U(:,:,j)'*Apr(:,j));
    end
    Apr2 = A(:,2:hedgesize) - A(:,1)*ones(1,hedgesize-1);
    [U(:,:,hedgesize),D,V]=svd(Apr2(:,[1:hedgesize-2]),'econ');
    b(:,hedgesize) = Apr2(:,hedgesize-1) - U(:,:,hedgesize)*(U(:,:,hedgesize)'*Apr2(:,j));
    
    segindex = 1:m;
    segindex(edgelist(i,1:hedgesize)) =[];
    C = data(:,segindex);
    Ctemp = bsxfun(@minus, C , A(:,hedgesize));
    for j = 1:hedgesize-1,
        Ctemp2 = Ctemp-U(:,:,j)*(U(:,:,j)'*Ctemp);
        dis(j,:) = sum(b(:,j).^2)-(b(:,j)'*Ctemp2).^2./sum(Ctemp2.^2);
        
    end
    Ctemp = bsxfun(@minus, C , A(:,1));
    Ctemp2 = Ctemp-U(:,:,hedgesize)*(U(:,:,hedgesize)'*Ctemp);
    dis(hedgesize,:) = sum(b(:,hedgesize).^2)-(b(:,hedgesize)'*Ctemp2).^2./sum(Ctemp2.^2);
    dis(hedgesize+1,:) = sum(Ctemp2.^2) -(b(:,hedgesize)'*Ctemp2).^2./sum(b(:,hedgesize).^2);
    directwedgelist(:,:,i) =dis;
end
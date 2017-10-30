function classes=directwedge_lap(directwedgelist,edgelist,c,N)
classes = zeros(1,N);
[edgenum,hyperedgesize]=size(edgelist);
A = zeros(N,N);
negw = 1/(hyperedgesize)/(hyperedgesize-1);
posw = 1/(hyperedgesize);

for i = 1:edgenum,
    segindex = 1:N;
    segindex(edgelist(i,:)) =[];
    A(edgelist(i,:),segindex) =  A(edgelist(i,:),segindex) + ones(hyperedgesize,1)*(posw*directwedgelist(hyperedgesize+1,:,i)...
        -negw*sum(directwedgelist(:,:,i)));
    for j = 1:hyperedgesize,
        A(edgelist(i,j),segindex) = A(edgelist(i,j),segindex)  + (posw+negw)*directwedgelist(j,:,i);
    end
end

A = A+A';
for i = 1:N,
    A(i,i) = 0;
end

A(A<0)=0;
D = diag(sum(A,2));

zerovlabel = find(diag(D)<=0);
classes(zerovlabel) = unidrnd(c,[1,length(zerovlabel)]);
nonzerovlabel = find(diag(D)>0);
A = A(nonzerovlabel,nonzerovlabel);
Dhalf  = D(nonzerovlabel,nonzerovlabel)^(-0.5);
NL = Dhalf*A*Dhalf ;
[U,D]=eigs(NL,c);
%U=U./(sqrt(sum(U.^2,2))*ones(1,c));
U = Dhalf*U;
[L,C] = kmeanspp(U',c);
%classes = L';
classes(nonzerovlabel) = L';
classes(zerovlabel) = unidrnd(c,[1,length(zerovlabel)]);
end
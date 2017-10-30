function [c1_DGC, c2_DGC] = InHCut(mutualinfoest,ksize)
[rows,n]=size(mutualinfoest);
Adj=zeros(n,n);
for i=2:n,
    for j=1:i-1,
        row = (i-1)*(i-2)/2+j;
        for k = 1:n,
            if k~=i && k~=j,
                Adj(i,j) = Adj(i,j) - mutualinfoest(row, k);
                Adj(i,k) = Adj(i,k) + mutualinfoest(row, k);
                Adj(j,k) = Adj(j,k) + mutualinfoest(row, k);
            end
        end
    end
end

Adj = Adj+ Adj';
Adj(Adj<0)=0;

D = diag(sum(Adj,2));
Dhalf =  D^(-1/2);
NL = eye(n) - Dhalf*Adj*Dhalf;

if ksize>=0,
    try
        [U,D]=eigs(NL,2,'sm');
        [sortedU,indexU]=sort(U(:,2));
    catch
        c1_DGC = randsample(n,ksize);
        c2_DGC = randsample(n,n-ksize);
        return;
    end
    
    tempcut1=sum(sum(Adj(indexU(1:ksize),indexU(ksize+1:n))));
    tempcut2=sum(sum(Adj(indexU(1:n-ksize),indexU(n-ksize+1:n))));
    if tempcut1<tempcut2,
        c1_DGC=indexU(1:ksize);
        c2_DGC=indexU(ksize+1:n);
    else
        c1_DGC=indexU(1:n-ksize);
        c2_DGC=indexU(n-ksize+1:n);
    end
    
else
    try
        [U,D]=eigs(NL,2,'sm');
        tempU = Dhalf*U;
        [sortedU,indexU]=sort(U(:,2));
    catch
        c1_DGC = randsample(n,n/2);
        c2_DGC = randsample(n,n/2);
        return;
    end
    
    ncutval=inf;
    iopt=0;
    for i=1:n-1,
        tempcut=sum(sum(Adj(indexU(1:i),indexU(i+1:n))));
        vol1=sum(sum(Adj(indexU(1:i),:)));
        vol2=sum(sum(Adj(indexU(i+1:n),:)));
        tempncutval=tempcut*(1/vol1+1/vol2);
        if ncutval>tempncutval,
            iopt=i;
            ncutval=tempncutval;
        end
    end
    
    c1_DGC=indexU(1:iopt);
    c2_DGC=indexU(iopt+1:n);
end
end
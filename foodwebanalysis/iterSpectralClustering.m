function [classes,objvalue,objcond,lap] = iterSpectralClustering(GM,c)
%spectral multiclass clustering via iteratively 2-partitioning
N = size(GM,1);
classes = zeros(1,N);
objvalue = zeros(c-1,1);
objcond = zeros(c-1,1);
lap = zeros(c-1,1);
effectlist = 1:N ;
for tempc = 1:c-1,
    D = diag(sum(GM,2));
    zerovlabel = find(diag(D)<=1e-12);
    classes(effectlist(zerovlabel)) = c+1;
    nonzerovlabel = find(diag(D)>1e-12);
    tempN= length(nonzerovlabel);
    GM = GM(nonzerovlabel,nonzerovlabel);
    D = D(nonzerovlabel,nonzerovlabel);
    NL = D^(-0.5)*GM*D^(-0.5);
 
    [U,D1]=eigs(NL,2);
    temp = U(:,1).*U(:,2);
    a = find(temp>0);
    b = find(temp<0);
    if U(a(1),1)*U(b(1),1)>0,
        h =2;
        lap(tempc) = 1-D1(2,2);
    else 
        h =1;
        lap(tempc) = 1-D1(1,1);
    end
    
    [sortedU,indexU]=sort(D^(-0.5)*U(:,h));
    ncutval=inf;
    iopt=0;
    for i=1:tempN-1,
        try
            tempcut=sum(sum(GM(indexU(1:i),indexU(i+1:tempN))));
        catch
            a = 1;
        end
        vol1=sum(sum(GM(indexU(1:i),:)));
        vol2=sum(sum(GM(indexU(i+1:length(nonzerovlabel)),:)));
        tempncutval=tempcut*(1/vol1+1/vol2);
        if ncutval>tempncutval,
            iopt=i;
            ncutval=tempncutval;
            objvalue(tempc) = ncutval;
            objcond(tempc) = tempcut*max(1/vol1,1/vol2);
        end
    end
    if iopt<=length(nonzerovlabel)/2,
        classes(effectlist(nonzerovlabel(indexU(1:iopt)))) = tempc;
        GM = GM(indexU(iopt+1:end),indexU(iopt+1:end));
        effectlist = effectlist(nonzerovlabel(indexU(iopt+1:end)));
    else
        classes(effectlist(nonzerovlabel(indexU(iopt+1:end)))) = tempc;
        GM = GM(indexU(1:iopt),indexU(1:iopt));
        effectlist = effectlist(nonzerovlabel(indexU(1:iopt)));
    end
end
classes(effectlist) = c;
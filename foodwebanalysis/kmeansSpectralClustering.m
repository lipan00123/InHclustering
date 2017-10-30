function classes = kmeansSpectralClustering(GM,c)
%spectral multiclass clustering via kmeans++
N = size(GM,1);
classes = zeros(1,N);
D = diag(sum(GM,2));
zerovlabel = find(diag(D)<=1e-12);
classes(zerovlabel) = unidrnd(c,[1,length(zerovlabel)]);
nonzerovlabel = find(diag(D)>1e-12);
GM = GM(nonzerovlabel,nonzerovlabel);
D = D(nonzerovlabel,nonzerovlabel);
NL = D^(-0.5)*GM*D^(-0.5); 
[U,D1]=eigs(NL,c);
U = U;
[L,C] = kmeanspp(U',c);
classes(nonzerovlabel) = L';
classes(zerovlabel) = unidrnd(c,[1,length(zerovlabel)]);
function [classes,errlabel,err_ols,timing]=directwedge_lapSS(data,c,p,edgenum,label,T)
%   Use inhomogeneous hypergraph to segment subspace
%
%   Input:
%
%   'data' - D*N data matrix where each column corresponds to a data point.
%            So there are N data points lying in a D-dim ambient space.
%
%   'c' - number of clusters
%
%   'p' - dim of the affine spaces
%
%   'edgenum' - number of pairs of vertices that will be sampled to compute
%               similarity. The sampling approach is similar to 
%               Spectral Curvature Clustering.
%
%   'label' - true labels, just used for computing error
%
%   'T' - number of iterations
%
%   Output:
%
%   'classes' - inferred labels of data points
%   
%   'errlabel' - number of errors in labels for each iteration
%
%   'err_ols' - distortion to the pure affine spaces via computing ols
%
%   'timing' - cpu time
%
%   If you have any questions please email lipan00123@gmail.com
%
%   Most Relevant Publications:
%   1. Inhomogeneous Hypergraph Clustering with Applications (NIPS2017)

[dim, N] = size(data);
templabel = ones(1, N);
errlabel = ones(T,1);
err_ols =  zeros(T,p+1);

minerr  =inf;
timing = zeros(T,1);
for t = 1:T,
    timing(t) = cputime;
    
    %sample pairs 
    edgelist=generate_hedgefordense(templabel,p,edgenum);
    
    %compute inhomogeneous distances 
    directwedgelist=generate_directwedgelist(data, edgelist);
    
    sortdirectwedgelist = max(directwedgelist,[],1);
    
    sortdirectwedgelist = sort(sortdirectwedgelist(:));
    
    %estimate sigma
    sigmalist = 10*sortdirectwedgelist(round(N*edgenum./c.^(1:p+1)));
    
    for sigmaindex = 1:length(sigmalist),
        
        %compute inhomogeneous weights
        tempdirectwedgelist = exp(-directwedgelist/sigmalist(sigmaindex));
        
        %inhomogeneous clustering
        classes=directwedge_lap(tempdirectwedgelist,edgelist,c,N);
        
        %compute ols to affine subspaces
        err_ols(t,sigmaindex) = errols(data,p,classes);
        if  minerr>err_ols(t,sigmaindex),
            minerr = err_ols(t,sigmaindex);
            templabel = classes;
        end
    end
    timing(t) = cputime-timing(t);
    
    %compute missclassfication errors 
    [errlabel(t),assignment] = missclassf(templabel,label);
end
timing=cumsum(timing);
err_ols = min(err_ols,[],2);
end
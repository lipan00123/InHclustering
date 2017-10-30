function mutualinfoest = samplelist(mutualinfoest,p);
    [pairs,n]= size(mutualinfoest);
    mask = zeros(pairs, n);
    for i = 3:n,
        for j = 2:i-1,
            for k = 1:j-1,
                if unifrnd(0,1)<p,
                    mask((i-1)*(i-2)/2+j,k) = 1;
                    mask((i-1)*(i-2)/2+k,j) = 1;
                    mask((j-1)*(j-2)/2+k,i) = 1;
                end
            end
        end
    end
    for i = 2:n,
        for j = 1:i-1,
            mask((i-1)*(i-2)/2+j,i) = 1;
            mask((i-1)*(i-2)/2+j,j) = 1;
        end
    end
    mutualinfoest = mutualinfoest.*mask;
end
function motiflist = checkmotif(Gmotif,method)
motiflist = [];
if strcmp(method,'D-bifan'),
    for i = 1:4,
        for j = 1:4,
            if j~=i,
                temp = 1:4;
                temp([i j]) = [];
                if sum(sum(Gmotif([i j],temp)))==4 && sum(sum(Gmotif(temp,[i j])))==0,
                    motiflist = [i j temp];
                    return;
                end
            end
        end
    end
end
end
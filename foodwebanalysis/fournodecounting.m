function [GM,mylist]= fournodecounting(G,Glist,motif,c)
N = length(Glist);
GM = zeros(N,N);
cur =1;
mylist = [];
Gundirlist=Glist;
Dundir =zeros(N,1);
if motif(1) == 'D',
    for i = 1:N,
        tempedge = Glist{i};
        for j = 1:length(tempedge),
            if G(tempedge(j),i) == 0,
                Gundirlist{tempedge(j)} = [ Gundirlist{tempedge(j)} i];
            end
        end
    end
    for i =1:N,
        Dundir(i) = length(Gundirlist{i});
    end
end

[sortedDundir, sortedindex] = sort(Dundir,'descend');
pos(sortedindex) = 1:N;

if strcmp(motif,'D-bifan')
    for i = 1:N,
        tempv = Gundirlist{sortedindex(i)};
        Gundirlist{sortedindex(i)} = [];
        for ini = 1:N,
            Ulist{ini} = [];
        end
        for j = 1:length(tempv),
            tempu = Gundirlist{tempv(j)};
            del = 0;
            for k = 1:length(tempu),
                if tempu(k)~=sortedindex(i),
                    Ulist{tempu(k)} = [Ulist{tempu(k)} tempv(j)];
                else
                    del = k;
                end
            end
            if del>0,
                Gundirlist{tempv(j)}(del) =[];
            end
        end
        for ini = 1:N,
            if length(Ulist{ini})>=2,
                %mylist{cur} = [sortedindex(i) ini Ulist{ini}];
                %cur = cur +1
                tempedge = Ulist{ini};
                for check1 = 2:length(tempedge),
                    for check2 = 1:check1-1,
                        templist = [sortedindex(i) ini tempedge(check1) tempedge(check2)];
                        motiflist = checkmotif(G(templist,templist),'D-bifan');
                        if length(motiflist)>0,
                            if G(templist(motiflist(1)),templist(motiflist(2)))...
                                    +G(templist(motiflist(2)),templist(motiflist(1))) == 1,
                                GM(templist(motiflist(1)),templist(motiflist(2))) = ...
                                    GM(templist(motiflist(1)),templist(motiflist(2))) + c(2);
                            else
                                GM(templist(motiflist(1)),templist(motiflist(2))) = ...
                                    GM(templist(motiflist(1)),templist(motiflist(2))) + c(1);
                            end
                            if G(templist(motiflist(3)),templist(motiflist(4)))...
                                    +G(templist(motiflist(4)),templist(motiflist(3))) == 1,
                                GM(templist(motiflist(3)),templist(motiflist(4))) = ...
                                    GM(templist(motiflist(3)),templist(motiflist(4))) + c(2);
                            else
                                GM(templist(motiflist(3)),templist(motiflist(4))) = ...
                                    GM(templist(motiflist(3)),templist(motiflist(4))) + c(1);
                            end                            
                            mylist =[mylist;templist(motiflist)];
                        end
                    end
                end
            end
        end
        i
    end
    GM = GM + GM';
end

if strcmp(motif,'bifan')
    for i = 1:N,
        tempv = Gundirlist{sortedindex(i)};
        Gundirlist{sortedindex(i)} = [];
        for ini = 1:N,
            Ulist{ini} = [];
        end
        for j = 1:length(tempv),
            tempu = Gundirlist{tempv(j)};
            del = 0;
            for k = 1:length(tempu),
                if tempu(k)~=sortedindex(i),
                    Ulist{tempu(k)} = [Ulist{tempu(k)} tempv(j)];
                else
                    del = k;
                end
            end
            if del>0,
                Gundirlist{tempv(j)}(del) =[];
            end
        end
        for ini = 1:N,
            if length(Ulist{ini})>=2,
                %mylist{cur} = [sortedindex(i) ini Ulist{ini}];
                %cur = cur +1
                tempedge = Ulist{ini};
                for check1 = 2:length(tempedge),
                    for check2 = 1:check1-1,
                        templist = [sortedindex(i) ini tempedge(check1) tempedge(check2)];
                        motiflist = checkmotif(G(templist,templist),'D-bifan');
                        if length(motiflist)>0,
                            GM(templist,templist) = GM(templist,templist)+1;
                        end
                    end
                end
            end
        end
        i
    end
    for i=1:N,
        GM(i,i)=0;
    end
end

end
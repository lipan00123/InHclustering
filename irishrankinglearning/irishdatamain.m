clear;
% load data
load meath.csv
[m,n]=size(meath);
n=n-1;
fullindex0 = [];
for i = 1:m,
    if isempty(find(meath(i,:)==-1,1)),
        fullindex0 = [fullindex0;i];
    end
end
m0 = length(fullindex0);

% randomly sample from the whole dataset
samplenum = 500;
times = 100;

totallvl=zeros(n,times);
totallvl2=zeros(n,times);

for count = 1:times,

    fullindex= fullindex0(randperm(m0,samplenum));
    m = samplenum;

    fullranking = meath(fullindex,1:n);
    fullranking = fullranking' + ones(n,m);

    fullrankingpos = zeros(n,m);
    for i = 1:m,
        fullrankingpos(fullranking(:,i),i)=(1:n)';
    end

    % estimate mutual info
    mutualinfoest = zeros(n*(n-1)/2,n);
    for i = 2:n,
        for j = 1:i-1,
            row = (i-1)*(i-2)/2+j;
            for k = 1:n,
                if k ~= i && k ~= j,
                    freqlist = zeros(2,n);
                    for p = 1:m,
                        freqlist(1+(fullrankingpos(i,p)>fullrankingpos(j,p)),fullrankingpos(k,p)) = freqlist(1+(fullrankingpos(i,p)>fullrankingpos(j,p)),fullrankingpos(k,p))+ 1;
                    end
                    freqlist = freqlist/m;
                    px = sum(freqlist,2);
                    py = sum(freqlist,1);
                    temp = freqlist.*log(freqlist./(px*py));
                    mutualinfoest(row,k) = sum(temp (~isnan(temp)));
                    if mutualinfoest(row,k)<0,
                        a=1;
                        mutualinfoest(row,k)=0;
                    end
                else
                    mutualinfoest(row,k) = 1;
                end
            end
        end
    end

        
    % 1st level InH 2-partition
    tempmutualinfoest = mutualinfoest;
    [c1_InH, c2_InH] = InHCut(tempmutualinfoest,3);
    if length(c1_InH)>length(c2_InH),
        temp = c1_InH;
        c1_InH = c2_InH;
        c2_InH = temp;
    end
    totallvl(c1_InH,count) = 1;
    remainlist = find(totallvl(:,count)==0);
    
    % 2nd level InH 2-partition
    temprowlist = generate_rowlist(remainlist);
    tempmutualinfoest = mutualinfoest(temprowlist,remainlist);
    [c1_InH, c2_InH] = InHCut(tempmutualinfoest,3);
    if length(c1_InH)>length(c2_InH),
        temp = c1_InH;
        c1_InH = c2_InH;
        c2_InH = temp;
    end
    totallvl(remainlist(c1_InH),count) = 2;
    remainlist = find(totallvl(:,count)==0);

    % 3rd level InH 2-partition
    temprowlist = generate_rowlist(remainlist);
    tempmutualinfoest = mutualinfoest(temprowlist,remainlist);
    [c1_InH, c2_InH] = InHCut(tempmutualinfoest,3);
    if length(c1_InH)>length(c2_InH),
        temp = c1_InH;
        c1_InH = c2_InH;
        c2_InH = temp;
    end
    totallvl(remainlist(c1_InH),count) = 3;
    remainlist = find(totallvl(:,count)==0);

    % 4th level InH 2-partition
    temprowlist = generate_rowlist(remainlist);
    tempmutualinfoest = mutualinfoest(temprowlist,remainlist);
    [c1_InH, c2_InH] = InHCut(tempmutualinfoest,3);
    if length(c1_InH)>length(c2_InH),
        temp = c1_InH;
        c1_InH = c2_InH;
        c2_InH = temp;
    end
    totallvl(remainlist(c1_InH),count) = 4;
    remainlist = find(totallvl(:,count)==0);

    % 1st level Apr 2-partition
    tempmutualinfoest = mutualinfoest;
    [c1_Apr, c2_Apr] = Apr(tempmutualinfoest,3,n-3);
    if length(c1_Apr)>length(c2_Apr),
        temp = c1_Apr;
        c1_Apr = c2_Apr;
        c2_Apr = temp;
    end
    totallvl2(c1_Apr,count) = 1;
    remainlist = find(totallvl2(:,count)==0);

    % 2nd level Apr 2-partition
    temprowlist = generate_rowlist(remainlist);
    tempmutualinfoest = mutualinfoest(temprowlist,remainlist);
    [c1_Apr, c2_Apr] = Apr(tempmutualinfoest,3,length(remainlist)-3);
    if length(c1_Apr)>length(c2_Apr),
        temp = c1_Apr;
        c1_Apr = c2_Apr;
        c2_Apr = temp;
    end
    totallvl2(remainlist(c1_Apr),count) = 2;
    remainlist = find(totallvl2(:,count)==0);

    % 3rd level Apr 2-partition
    temprowlist = generate_rowlist(remainlist);
    tempmutualinfoest = mutualinfoest(temprowlist,remainlist);
    [c1_Apr, c2_Apr] = Apr(tempmutualinfoest,3,length(remainlist)-3);
    if length(c1_Apr)>length(c2_Apr),
        temp = c1_Apr;
        c1_Apr = c2_Apr;
        c2_Apr = temp;
    end
    totallvl2(remainlist(c1_Apr),count) = 3;
    remainlist = find(totallvl2(:,count)==0);

    % 4th level Apr 2-partition
    temprowlist = generate_rowlist(remainlist);
    tempmutualinfoest = mutualinfoest(temprowlist,remainlist);
    [c1_Apr, c2_Apr] = Apr(tempmutualinfoest,3,length(remainlist)-3);
    if length(c1_Apr)>length(c2_Apr),
        temp = c1_Apr;
        c1_Apr = c2_Apr;
        c2_Apr = temp;
    end
    totallvl2(remainlist(c1_Apr),count) = 4;
    remainlist = find(totallvl2(:,count)==0);
end

save('irish','totallvl','totallvl2')

%% compute successful rate
[n,times,mnum]=size(totallvl);
successInH = zeros(4,mnum);
successApr = zeros(4,mnum);
set1true = [1;4;13];
set2true = [2;5;6];
set3true = [7;8;9];

for i=1:times,
    for rindex = 1:mnum,
    set0 = find(totallvl(:,i,rindex)==0);
    set1 = find(totallvl(:,i,rindex)==1);
    set2 = find(totallvl(:,i,rindex)==2);
    set3 = find(totallvl(:,i,rindex)==3);
    set4 = find(totallvl(:,i,rindex)==4);
    temp1 = (sum((set0==set1true))==3)+(sum((set1==set1true))==3)+(sum((set2==set1true))==3)+(sum((set3==set1true))==3);
    temp2 = (sum((set0==set2true))==3)+(sum((set1==set2true))==3)+(sum((set2==set2true))==3)+(sum((set3==set2true))==3);
    temp3 = (sum((set0==set3true))==3)+(sum((set1==set3true))==3)+(sum((set2==set3true))==3)+(sum((set3==set3true))==3);
    
    successInH(1,rindex)=successInH(1,rindex)+temp1;
    successInH(2,rindex)=successInH(2,rindex)+temp2;
    successInH(3,rindex)=successInH(3,rindex)+temp3;
    successInH(4,rindex)=successInH(4,rindex)+(temp1+temp2+temp3==3);
    
    set0 = find(totallvl2(:,i,rindex)==0);
    set1 = find(totallvl2(:,i,rindex)==1);
    set2 = find(totallvl2(:,i,rindex)==2);
    set3 = find(totallvl2(:,i,rindex)==3);
    set4 = find(totallvl2(:,i,rindex)==4);
    temp1 = (sum((set0==set1true))==3)+(sum((set1==set1true))==3)+(sum((set2==set1true))==3)+(sum((set3==set1true))==3);
    temp2 = (sum((set0==set2true))==3)+(sum((set1==set2true))==3)+(sum((set2==set2true))==3)+(sum((set3==set2true))==3);
    temp3 = (sum((set0==set3true))==3)+(sum((set1==set3true))==3)+(sum((set2==set3true))==3)+(sum((set3==set3true))==3);
    
    successApr(1,rindex)=successApr(1,rindex)+temp1;
    successApr(2,rindex)=successApr(2,rindex)+temp2;
    successApr(3,rindex)=successApr(3,rindex)+temp3;
    successApr(4,rindex)=successApr(4,rindex)+(temp1+temp2+temp3==3);
    end
end

successInH = successInH/100;
successApr = successApr/100;

fprintf('The number of samples is %d\n', samplenum);
fprintf('The success rates of InH are: F.F. %.2f, F.G. %.2f, Ind. %.2f, overall %.2f\n', successInH(1), successInH(2), successInH(3), successInH(4));
fprintf('The success rates of Apr are: F.F. %.2f, F.G. %.2f, Ind. %.2f, overall %.2f\n', successApr(1), successApr(2), successApr(3), successApr(4));



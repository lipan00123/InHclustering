clear;

%extract data
str='1R2RC';
tempstr = strcat('./Hopkins155/',str,'/',str,'_truth');
load(tempstr);
class1 = find(s==1);
data1 = zeros(2*frames,length(class1));
for i = 1:length(class1),
    for j = 1:frames,
        data1(2*j-1:2*j,i)=x(1:2,class1(i),j);
    end
end
class2 = find(s==2);
data2 = zeros(2*frames,length(class2));
for i = 1:length(class2),
    for j = 1:frames,
        data2(2*j-1:2*j,i)=x(1:2,class2(i),j);
    end
end
class3 = find(s==3);
data3 = zeros(2*frames,length(class3));
for i = 1:length(class3),
    for j = 1:frames,
        data3(2*j-1:2*j,i)=x(1:2,class3(i),j);
    end
end
data = [data1 data2 data3];
label = [ones(1,length(class1)) 2*ones(1,length(class2)) 3*ones(1,length(class3))];
n = length(label);

%hyperedge size
hedgesize =5;
%subspace dim
p = 4;
%the number of pairs
edgenum = 300;
%The number of iterations to resample the pairs 
T =10;


classes=directwedge_lapSS(data,3,p,edgenum,label,T);

fprintf('data: %s; #points: %d, #sample pairs: %d\n', str, n, edgenum);
for i=1:T,
    fprintf('The number of missclassified points is %d in ith iteration\n', classes(i));
end




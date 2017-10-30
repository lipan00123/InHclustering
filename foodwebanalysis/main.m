clear;

%% read data
load Florida-bay.txt
dep = readtable('./Florida-bay-meta.csv');
dep = table2cell(dep);
m = size(Florida_bay,1);
V = max(max(Florida_bay))+1;

%% construct projected graph
% construct the Ajacency matrix %
G=zeros(V,V);
for i = 1:m,
    G(Florida_bay(i,1)+1,Florida_bay(i,2)+1) = 1;
end
G = G(1:V-1,1:V-1);

for i = 1:V-1,
    newGlist{i}=[];
    for j = 1:V-1,
        if G(i,j)==1,
            newGlist{i} =[newGlist{i} j];
        end
    end
end

%count the bifan motifs and make projection: GM is the projected graph %
[GM,mylist]= fournodecounting(G,newGlist,'D-bifan',[1 1]);
save(strcat('bifanno128_formal'),'GM','mylist');

%% spectral clustering
load('bifanno128_formal');
% first level
classes = kmeansSpectralClustering(GM,3);

% identify singleton in 1st level
sig_l1 =find(classes==4);
% identify three clusters
class1 = find(classes==classes(1)); 
class2 = find(classes==classes(9));
class3 = find(classes==classes(99));

% second level
% further partition class2
G2 = G(class2,class2);
for i = 1:length(class2),
    newG2list{i}=[];
end

for i = 1:length(class2),
    for j = 1:length(class2),
        if G2(i,j)==1,
            newG2list{i} =[newG2list{i} j];
        end
    end
end
%count the bifan motifs and make projection: GM2 is the projected subgraph of class2%
GM2= fournodecounting(G2,newG2list,'D-bifan',[1 1]);
classes2dep = iterSpectralClustering(GM2,2);
% identify singleton of class2 of 2nd level
sig_l2 =find(classes2dep==3);

% third level
% further partition class3
G3 = G(class3,class3);
for i = 1:length(class3),
    newG3list{i}=[];
end

for i = 1:length(class3),
    for j = 1:length(class3),
        if G3(i,j)==1,
            newG3list{i} =[newG3list{i} j];
        end
    end
end
%count the bifan motifs and make projection: GM2 is the projected subgraph of class3%
GM3= fournodecounting(G3,newG3list,'D-bifan',[1 1]);
classes3dep = iterSpectralClustering(GM3,2);
% identify singleton of class3 of 2nd level
sig_l3 =find(classes3dep==3);

%store the clustering results
save('foodwebfinalanalysis','class1','class2','class3','classes2dep','classes3dep','sig_l1','sig_l2','sig_l3');

%% draw figures
load foodwebfinalanalysis

classsort = [class1 class2(find(classes2dep==2)) class2(find(classes2dep==1)) class3(find(classes3dep==2)) class3(find(classes3dep==1))];
len(1) = length(class1);
len(2) = length(class2(find(classes2dep==2)));
len(3) = length(class2(find(classes2dep==1)));
len(4) = length(class3(find(classes3dep==2)));
len(5) = length(class3(find(classes3dep==1)));
classlabel = [ones(1,len(1)) 2*ones(1,len(2)) 3*ones(1,len(3)) 4*ones(1,len(4)) 5*ones(1,len(5))];


newclass2 = class2(find(classes2dep==2));
newclass3 = class2(find(classes2dep==1));
newclass4 = class3(find(classes3dep==2));
newclass5 = class3(find(classes3dep==1));

DOC(1)=find(classsort==123);
DOC(2)=find(classsort==124);
DOC(3)=find(classsort==125);
%DOC(4)=find(classsort==126);
DOC(4)=find(classsort==127);
classsort(DOC) = [];
classlabel(DOC) = [];

Gfigure = digraph(G(class1,class1));
figure(1);
h=plot(Gfigure,'NodeLabel',dep(class1,2),'Layout','layered');
XData1 = h.XData;
YData1 = h.YData;

Gfigure = digraph(G(newclass2,newclass2));
%figure(2);
h=plot(Gfigure,'NodeLabel',dep(newclass2,2),'Layout','layered');
XData2 = h.XData;
YData2 = h.YData;

Gfigure = digraph(G(newclass3,newclass3));
%figure(3);
h=plot(Gfigure,'NodeLabel',dep(newclass3,2),'Layout','layered');
XData3 = h.XData;
YData3 = h.YData;

Gfigure = digraph(G(newclass4,newclass4));
%figure(4);
h=plot(Gfigure,'NodeLabel',dep(newclass4,2),'Layout','layered');
XData4 = h.XData;
YData4 = h.YData;

Gfigure = digraph(G(newclass5,newclass5));
%figure(5);
h=plot(Gfigure,'NodeLabel',dep(newclass5,2),'Layout','layered');
XData5 = h.XData;
YData5 = h.YData;
basicnodecolor = h.NodeColor;

gap = 20;
myXData = [XData1/max(XData1)*80 XData2/max(XData2)*80 XData3/max(XData3)*80 XData4/max(XData4)*80 XData5/max(XData5)*80];
myYData = [YData1 YData2+gap YData3+2*gap YData4+3*gap YData5+4*gap];
mynodecolor = [ones(len(1),1)*[0 0.5 0];ones(len(2),1)*[0.8 0.8 0];ones(len(3),1)*[1 0 0];...
    ones(len(4),1)*[0 0 1];ones(len(5),1)*[0.8 0 0.8];];

mynodecolor(DOC,:) = [];
myXData(DOC) =[];
myYData(DOC) =[];

Gfigure = digraph(G(classsort,classsort));
h=plot(Gfigure,'NodeLabel',dep(classsort,2),'XData',myXData,'YData',myYData);
basicedgecolor = h.EdgeColor;
basiclinewidth = h.LineWidth;
edgelist = Gfigure.Edges.EndNodes;
basicnodecolor = h.NodeColor;

inversepair = [];
inversepos = [];
for i = 1:length(classsort)-1,
    for j = i+1:length(classsort),
        if G(classsort(i),classsort(j)) == 0 && G(classsort(j),classsort(i)) == 1 && classlabel(j)~=classlabel(i),
            inversepair= [inversepair;i j];
            temp = find(edgelist(:,2)==i);
            
            inversepos = [inversepos;temp(find(edgelist(temp,1)==j))];
        end
    end
end

myedgecolor=ones(length(edgelist(:,1)),1)*basicedgecolor;
myedgecolor(inversepos,:) = ones(length(inversepos),1)*[0 0 0];
mylinewidth=ones(length(edgelist(:,1)),1)*basiclinewidth;
mylinewidth(inversepos,:)=2;
figure(1);
h=plot(Gfigure,'MarkerSize',5,'NodeColor',mynodecolor,'XData',myYData,'YData',-myXData,'EdgeColor',myedgecolor,'LineWidth',mylinewidth);


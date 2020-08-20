function [convergence_max,convergence_min,q] = Containmentnode(adj)
%UNTITLED 此处显示有关此函数的摘要
%   此处显示详细说明


[q,e,r]=extract_connected_component(adj);
% l=BFSDivideGraph(adj);
a=0.2;
b=0.8;
for i=1:size(e,2)
p(i).c=ClusteringCoefficient(q(i).matrix);%子图的聚集数
p(i).d=diag(DU(q(i).matrix))';%子图的度矩阵
p(i).dc=a*mapminmax(p(i).c,0,1)+ b*mapminmax(p(i).d,0,1);               %子图度和聚集数组合评价 归一化求加权
[~,x(i)]  =max(p(i).dc);%最大聚集数及位于子图中的位置
[~,y(i)]  =min(p(i).dc);%最小聚集数及位于子图中的位置
convergence_max(i)=q(i).kept(x(i)); %子图中具有最大聚集数的个体
convergence_min(i)=q(i).kept(y(i)); %子图中具有最小聚集数的个体

end
convergence_max;
convergence_min;
end


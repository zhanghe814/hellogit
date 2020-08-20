function C=ClusteringCoefficient(a)
n=size(a,2);%a为邻接矩阵
C=zeros(1,n);
for i=1:n
   aa=find(a(i,:)==1);%与节点i相连的节点集合
   if isempty(aa)
       C(i)=0;
   else
       m=length(aa);%i的度或者临域维数
       if m==1
           C(i)=0;
       else
           B=a(aa,aa);%B是i的邻域
           C(i)=length(find(B==1))/m/(m-1);%B==0时表示与i连接的两点之间无连接
       end
   end
end
C;
cc=mean(C);
    
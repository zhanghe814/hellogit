function C=ClusteringCoefficient(a)
n=size(a,2);%aΪ�ڽӾ���
C=zeros(1,n);
for i=1:n
   aa=find(a(i,:)==1);%��ڵ�i�����Ľڵ㼯��
   if isempty(aa)
       C(i)=0;
   else
       m=length(aa);%i�ĶȻ�������ά��
       if m==1
           C(i)=0;
       else
           B=a(aa,aa);%B��i������
           C(i)=length(find(B==1))/m/(m-1);%B==0ʱ��ʾ��i���ӵ�����֮��������
       end
   end
end
C;
cc=mean(C);
    
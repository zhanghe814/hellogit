function [convergence_max,convergence_min,q] = Containmentnode(adj)
%UNTITLED �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��


[q,e,r]=extract_connected_component(adj);
% l=BFSDivideGraph(adj);
a=0.2;
b=0.8;
for i=1:size(e,2)
p(i).c=ClusteringCoefficient(q(i).matrix);%��ͼ�ľۼ���
p(i).d=diag(DU(q(i).matrix))';%��ͼ�ĶȾ���
p(i).dc=a*mapminmax(p(i).c,0,1)+ b*mapminmax(p(i).d,0,1);               %��ͼ�Ⱥ;ۼ���������� ��һ�����Ȩ
[~,x(i)]  =max(p(i).dc);%���ۼ�����λ����ͼ�е�λ��
[~,y(i)]  =min(p(i).dc);%��С�ۼ�����λ����ͼ�е�λ��
convergence_max(i)=q(i).kept(x(i)); %��ͼ�о������ۼ����ĸ���
convergence_min(i)=q(i).kept(y(i)); %��ͼ�о�����С�ۼ����ĸ���

end
convergence_max;
convergence_min;
end


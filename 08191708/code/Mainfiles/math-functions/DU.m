function D = DU(A)
%UNTITLED2 �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
[n,~]=size(A);
for i=1:n
    for j=1:n
        m(i)=sum(A(i,:));
    end
end

D=diag(m);
end


function D = DU(A)
%UNTITLED2 此处显示有关此函数的摘要
%   此处显示详细说明
[n,~]=size(A);
for i=1:n
    for j=1:n
        m(i)=sum(A(i,:));
    end
end

D=diag(m);
end


clc;
clear all;
close all;
% w_matrix=load('w_matrix.txt');
% w_matrix=w_matrix((1:end),(1:30));
qianzhi=textread('qianzhi.txt');
qianzhi_log=(qianzhi~=0);

for i=1:size(qianzhi,1)
    n(i)=sum(qianzhi_log(i,:));
end
plot(n)
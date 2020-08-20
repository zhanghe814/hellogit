clc;
clear all;
close all;
% w_matrix=load('w_matrix.txt');
% w_matrix=w_matrix((1:end),(1:30));
% qianzhi=textread('qianzhi.txt');
% qianzhi_log=(qianzhi~=0);

% for i=1:size(qianzhi,1)
%     n(i)=sum(qianzhi_log(i,:));
% end
% plot(n)
load('data2020-08-19_20-07-23.mat')
size(data,1);
for j=1:size(data(1).agents{1}(:,1),1)
   for i=1:size(data,1)
       follow.a{j}(i,:)=data(i).agents{1}(j,:);
   end
end

   for i=1:size(data,1)
      leader.a(i,:)=data(i).gamma(1,:);
   end
% %%%计算邻接矩阵
% for k=1:1:size(data(1).agents{1}(:,1),1)
%     for j=1:1:size(data(1).agents{1}(:,1),1)
%         
%           distanceX{k}=follow.a{k}(:,1)-follow.a{j}(:,1);
%           distanceY{k}=follow.a{k}(:,2)-follow.a{j}(:,2);
%     end
% end
% r=1.5*5;
% epsilon=0.1;
% for k=size(follow.a{1}(:,1),1)
% for i=1:1:size(data(1).agents{1}(:,1),1)
%     for j=1:1:size(data(1).agents{1}(:,1),1)
%         
%         w_matrix(i,j)=adj_matrix([distanceX{i}(j,:),distanceY{i}(j,:)],r);
%     end
% end
% 
% end

figure(1)%x轴位置
for k=1:1:size(data(1).agents{1}(:,1),1)
    plot(follow.a{k}(:,1))
    hold on
end
plot(leader.a(1,:))

figure(2)%y轴位置
for k=1:1:size(data(1).agents{1}(:,1),1)
    plot(follow.a{k}(:,2))
    hold on
end
plot(leader.a(2,:))

figure(3)%x轴速度
for k=1:1:size(data(1).agents{1}(:,1),1)
    plot(follow.a{k}(:,3))
    hold on
    plot(leader.a(3,:),'o')
end
plot(leader.a(3,:))
figure(4)%y轴速度
for k=1:1:size(data(1).agents{1}(:,1),1)
    plot(follow.a{k}(:,4))
    hold on
    plot(leader.a(4,:),'o')
end
plot(leader.a(4,:))
figure(5)%xy轴位置
for k=1:1:size(data(1).agents{1}(:,1),1)
    plot(follow.a{k}(:,1),follow.a{k}(:,2))
    hold on
end
plot(leader.a(1,:),leader.a(2,:))

% hold on
% r=10; theta=0:pi/100:2*pi;
% x=r*cos(theta); y=r*sin(theta);
% plot(x,y,'-')
% fill(x,y,'c')
% hold on 
% r=4;
% x=r*cos(theta)+30; y=r*sin(theta)+3;
% plot(x,y,'-')
% fill(x,y,'c')
% hold on
% r=5;
% x=r*cos(theta)+100; y=r*sin(theta);
% plot(x,y,'-')
% fill(x,y,'c')



function out = odeRhs(~,odeVec,settings,simObj) 
%ODERHS returns right hand side of 1st order ODE "dp/dt = u(t,p)"
%   where in our case u_i(t,p)is as on page 9 of the paper
% global qianyin  count ;
% qianyin={};
%Copy settings into local variables
NSpecies = settings.NSpecies;
NAgents = settings.NAgents;
TotAgent = size(odeVec,1) / 4;
agents = reshape(odeVec,TotAgent,4);
r = settings.r;
d = settings.d;
Species_c1 = settings.Species_c1;
Species_c2 = settings.Species_c2;
mass = settings.mass;
w_matrix=zeros(TotAgent);
%range_lower_limit = -10;
%range_upper_limit = 10;


[mesh1X, mesh2X] = meshgrid(agents(:,1),agents(:,1));
[mesh1Y, mesh2Y] = meshgrid(agents(:,2),agents(:,2));
distanceMatrixX = mesh2X-mesh1X; % distance matrix in x direction
distanceMatrixY = mesh2Y-mesh1Y; % distance matirx in y direction
distanceCellX = mat2cell(distanceMatrixX,NAgents,NAgents);
distanceCellY = mat2cell(distanceMatrixY,NAgents,NAgents);

[mesh1VX, mesh2VX] = meshgrid(agents(:,3),agents(:,3));
[mesh1VY, mesh2VY] = meshgrid(agents(:,4),agents(:,4));
velocityMatrixX = mesh2VX-mesh1VX; % distance matrix in x direction
velocityMatrixY = mesh2VY-mesh1VY; % distance matirx in y direction
velocityCellX = mat2cell(velocityMatrixX,NAgents,NAgents);
velocityCellY = mat2cell(velocityMatrixY,NAgents,NAgents);

%-------------------------------Calc u_alpha-----------------------------%
u_alpha_temp = cell(size(distanceCellX));
u_alpha_temp_2 = cell(length(distanceCellX),1);
   
%  %%%%%计算邻接矩阵
%  for i = 1:length(distanceCellX{1,1}(1,:))
% %             fi_g = zeros(1,2); %Gradient based term
% %             fi_d = zeros(1,2); %velocity consensus term
%             for j = 1:length(distanceCellX{1,1}(:,1))
%                % if k == l    
%                     if ((j~=i) && (norm([distanceCellX{1,1}(j,i),distanceCellY{1,1}(j,i)]) < r(1)))
%                      w_matrix(i,j)=adj_matrix([distanceCellX{1,1}(j,i),distanceCellY{1,1}(j,i)],r(1));
%                     end
%                 
%             end
%  end
%     w_matrix(w_matrix~=0)=1;
%     [convergence_max,q]=Containmentnode(w_matrix);
    
    
    
% flag_1=1;
    for k = 1:length(distanceCellX)
    for l = 1:length(distanceCellX)
        for i = 1:length(distanceCellX{l,k}(1,:))
            fi_g = zeros(1,2); %Gradient based term
            fi_d = zeros(1,2); %velocity consensus term
            for j = 1:length(distanceCellX{l,k}(:,1))
                if k == l    
                    if ((j~=i) && (norm([distanceCellX{l,k}(j,i),distanceCellY{l,k}(j,i)]) < r(k)))
%                         for n=1:size(convergence_max,2)
%                             if j==convergence_max(n)
%                                 flag_1=0;
%                             end
%                             
%                         end
%                         if flag_1==0
%                             fi_d(1,:)=[0,0];
%                         else
                        fi_g(1,:) = fi_g + phi_alpha(sigma_norm([distanceCellX{l,k}(j,i),distanceCellY{l,k}(j,i)]),d(k),r(k))*gradient_sigma_norm([distanceCellX{l,k}(j,i),distanceCellY{l,k}(j,i)]);
                        fi_d(1,:) = fi_d + adj_matrix([distanceCellX{l,k}(j,i),distanceCellY{l,k}(j,i)],r(k))*[velocityCellX{l,k}(j,i),velocityCellY{l,k}(j,i)];
                        w_matrix(i,j)=adj_matrix([distanceCellX{l,k}(j,i),distanceCellY{l,k}(j,i)],r(k));
%                         end
                    end
                else
                    if  (norm([distanceCellX{l,k}(j,i),distanceCellY{l,k}(j,i)]) < r(k))
                        fi_g(1,:) = fi_g + phi_alpha(sigma_norm([distanceCellX{l,k}(j,i),distanceCellY{l,k}(j,i)]),d(k),r(k))*gradient_sigma_norm([distanceCellX{l,k}(j,i),distanceCellY{l,k}(j,i)]);
                        fi_d(1,:) = fi_d + adj_matrix([distanceCellX{l,k}(j,i),distanceCellY{l,k}(j,i)],r(k))*[velocityCellX{l,k}(j,i),velocityCellY{l,k}(j,i)];
                        w_matrix(i,j)=adj_matrix([distanceCellX{l,k}(j,i),distanceCellY{l,k}(j,i)],r(k));
                    end
                end
            end
            fi_g = Species_c1(l,k)*fi_g;
            fi_d = Species_c2(l,k)*fi_d;
            u_alpha_temp{l,k}(i,:) = fi_g + fi_d;
        end
    end
    end



for k = 1:length(u_alpha_temp(1,:))
    temp = zeros(size(u_alpha_temp{1,k}));
    for j = 1:length(u_alpha_temp(:,1))
        temp =  temp + u_alpha_temp{j,k};
    end
    u_alpha_temp_2{k,1} = temp;
end

u_alpha = cell2mat(u_alpha_temp_2);
  
    
if strcmp(settings.SimMode, 'algorithm1')
    u = u_alpha;

elseif strcmp(settings.SimMode, 'algorithm2') || strcmp(settings.SimMode, 'algorithm3')
    
    %copy settings into local variables
    gamma_agent = simObj.gamma;
    qd = gamma_agent(1:2);
    pd = gamma_agent(3:4);
    Gamma_c1 = settings.Gamma_c1;
    Gamma_c2 = settings.Gamma_c2;
    random = settings.randomBool;
    force_mag = 3;
    
% % % ----------------------------Calc u_gamma------------------------------%    
    
%     %Agent pos and vel -  gamma agent pos and vel
%     delta_pos = [(agents(:,1) - repmat(qd(1),TotAgent,1)),(agents(:,2) - repmat(qd(2),TotAgent,1))];
%     delta_vel = [(agents(:,3) - repmat(pd(1),TotAgent,1)),(agents(:,4) - repmat(pd(2),TotAgent,1))];
%    
%     %Repeat gamma constants (NAgents) number of times
%     Gamma_c1_rep = repmat(repeat(NAgents,Gamma_c1)',1,2);
%     Gamma_c2_rep = repmat(repeat(NAgents,Gamma_c2)',1,2);
%     
%     u_gamma = -Gamma_c1_rep.*delta_pos - Gamma_c2_rep.*delta_vel;
     
  

%%%%%%%%%%改为牵引运动
    w_matrix(w_matrix~=0)=1;
%     fid=fopen('w_matrix.txt','a');
%     fprintf(fid,'%d ',w_matrix);
%     fprintf(fid,'\r\n');
%     fclose(fid);
%     
%     save('w_matrix.mat','w_matrix') ;
  
    [convergence_max,convergence_min,q]=Containmentnode(w_matrix);
    
    fid=fopen('qianzhi.txt','a');
    fprintf(fid,'%d ',convergence_max);
    fprintf(fid,'\r\n');
    fclose(fid);
    
%     cc=ClusteringCoefficient_mean(w_matrix);
    for i=1:size(convergence_max,2) 
        %牵制节点运动 Agent pos and vel（牵制节点） -  gamma agent pos and vel（虚拟领导者）
        delta_pos(convergence_max(i),:) = [(agents(convergence_max(i),1) - repmat(qd(1),1,1)),(agents(convergence_max(i),2) - repmat(qd(2),1,1))];
        delta_vel(convergence_max(i),:) = [(agents(convergence_max(i),3) - repmat(pd(1),1,1)),(agents(convergence_max(i),4) - repmat(pd(2),1,1))];
        
    %%%%非牵制节点的
    non_qianzhi_pos=find(q(i).kept~=convergence_max(i));
   
    for k=1:size(non_qianzhi_pos,1)
    non_qianzhi(k)=q(i).kept(non_qianzhi_pos(k));
    delta_pos(non_qianzhi(k),:) =0; %[(agents(non_qianzhi(k),1) - agents(convergence_max(i),1)),(agents(non_qianzhi(k),2) - agents(convergence_max(i),2))];%[agents(convergence_max(i),1),agents(convergence_max(i),2)];%[(agents(non_qianzhi(k),1) - agents(convergence_max(i),1)),(agents(non_qianzhi(k),2) - agents(convergence_max(i),2))];
    delta_vel(non_qianzhi(k),:) =0; %[(agents(non_qianzhi(k),3) - agents(convergence_max(i),3)),(agents(non_qianzhi(k),4) - agents(convergence_max(i),4))];%[agents(convergence_max(i),3),agents(convergence_max(i),4)];%[(agents(non_qianzhi(k),3) - agents(convergence_max(i),3)),(agents(non_qianzhi(k),4) - agents(convergence_max(i),4))];
    end
        
    end
    
%     cc=ClusteringCoefficient_mean(w_matrix);
%     D=diag(DU(w_matrix));
%     D_mean=mean(D);
%     ccc=cc+D_mean;
%     fid=fopen('ccc.txt','a');
%     fprintf(fid,'%d ',ccc);
%     fprintf(fid,'\r\n');
%     fclose(fid);
    
    %%%%%%解决聚集度小问题  求度和聚集度加权最小的也获得目标信息
    if strcmp(settings.SimMode, 'algorithm3')
    cc=ClusteringCoefficient_mean(w_matrix);
    D=diag(DU(w_matrix));
    D_mean=mean(D);
    cc=cc+D_mean;
    fid=fopen('cc.txt','a');
    fprintf(fid,'%d ',cc);
    fprintf(fid,'\r\n');
    fclose(fid);
    if size(convergence_max,2) <3
        if cc<6%6
        for i=1:size(convergence_min,2)
           delta_pos(convergence_min(i),:) = [(agents(convergence_min(i),1) - repmat(qd(1),1,1)),(agents(convergence_min(i),2) - repmat(qd(2),1,1))];
           delta_vel(convergence_min(i),:) = [(agents(convergence_min(i),3) - repmat(pd(1),1,1)),(agents(convergence_min(i),4) - repmat(pd(2),1,1))];
        end
        
        end
            
    end
    end
    
    
%      if (size(non_qianzhi_pos,1))>=(size(w_matrix,1)-1)
%         
%         delta_pos= [(agents(:,1) - repmat(qd(1),TotAgent,1)),(agents(:,2) - repmat(qd(2),TotAgent,1))];
%         delta_vel= [(agents(:,3) - repmat(pd(1),TotAgent,1)),(agents(:,4) - repmat(pd(2),TotAgent,1))];   
%     end
        %Repeat gamma constants (NAgents) number of times
%     if (size(convergence_max,2)<=3)
%         Gamma_c1=4*Gamma_c1;
%         Gamma_c2=4*Gamma_c2;
%     end

    Gamma_c1_rep = repmat(repeat(NAgents,Gamma_c1)',1,2);
    Gamma_c2_rep = repmat(repeat(NAgents,Gamma_c2)',1,2);
    
%     u_gamma = -Gamma_c1_rep.*delta_pos - Gamma_c2_rep.*delta_vel;
    
    
%     u_gamma = -TotAgent/size(convergence_max,2)*Gamma_c1_rep.*delta_pos - TotAgent/size(convergence_max,2)*Gamma_c2_rep.*delta_vel;


   %%%调节稳定性
     gamma_i=TotAgent/size(convergence_max,2);
     gamma_d=TotAgent/size(convergence_max,2);
     if gamma_i>0.8*TotAgent
         gamma_i=0.8*TotAgent;
     end
     if gamma_d>0.8*TotAgent
         gamma_d=0.8*TotAgent;
     end
     u_gamma = -2.8*gamma_i*Gamma_c1_rep.*delta_pos - 0.82*gamma_d*Gamma_c2_rep.*delta_vel;
  
    
    
    
    
 
%-------------------------Calc u_rand--------------------------------%

if random
    distanceCellXY = mat2cell(sqrt(distanceMatrixX.^2 + distanceMatrixY.^2),NAgents,NAgents);
    for i = 1:length(distanceCellXY)
        for j = 1:length(distanceCellXY)
            distanceCellXY{j,i} = distanceCellXY{j,i} < r(i);
        end
    end
    
    neighbors_count = repmat((sum(cell2mat(distanceCellXY)~=0,1)' > 7),1,2); %Agents with more than 6 neighbors (7 since 1 itself) repeated columnwise for x and y forces
    u_rand = zeros(size(agents(:,1:2)));
    for i = 1:length(u_rand(:,1))
        if neighbors_count(i,:) ~= 0
            u_rand(i,:) = force_mag*(agents(i,1:2)-qd)/norm(agents(i,1:2)-qd);
        end
    end
else
    u_rand = 0;
end

u = u_alpha + u_gamma + u_rand;

else
    error('Set mode to either algorithm1, algorithm2 or algorithm3');
end
    
if strcmp(settings.SimMode, 'algorithm3')
  %Calculate and add u_beta if algorithm 3
  obstacles = simObj.obstacles;
  NObstacles = length(obstacles(:,1));
  r_prime = settings.r_prime;
  d_prime = settings.d_prime;
  Beta_c1 = settings.Beta_c1;
  Beta_c2 = settings.Beta_c2;
  Agents = simObj.agents; %Different from agents used above. Above is array.
  u_beta = cell(size(Agents));
  
  D=diag(DU(w_matrix));
  for i = 1:length(Agents)
      for j = 1:length(Agents{1,i}(:,1))
          fi_g = zeros(1,2); %Gradient based term
          fi_d = zeros(1,2); %velocity consensus term
          for k = 1:NObstacles
              if norm(Agents{1,i}(j,1:2) - obstacles(k,1:2)) <= (r_prime(i)+obstacles(k,3)) %vicinity check. If obstacle within agent's obstacle sensing range
                  mu_beta = obstacles(k,3)/norm(Agents{1,i}(j,1:2) - obstacles(k,1:2)); 
                  a_beta = (Agents{1,i}(j,1:2) - obstacles(k,1:2))/norm(Agents{1,i}(j,1:2) - obstacles(k,1:2));
                  P_beta = eye(2) - a_beta'*a_beta;
                  q_beta = mu_beta*Agents{1,i}(j,1:2) + (1-mu_beta)*obstacles(k,1:2);
                  p_beta = (mu_beta*P_beta*Agents{1,i}(j,3:4)')';
                 %%%%%解决避障问题 
                  if (D(i)<6)
                      D(i)=5;
                  end
                  if(D(i)>0.5*30)
                      D(i)=0.5*30; 
                  end
%                   fi_g = fi_g + Beta_c1(i)*phi_alpha(sigma_norm(q_beta - Agents{1,i}(j,1:2)),d_prime(i),r_prime(i))*gradient_sigma_norm(q_beta - Agents{1,i}(j,1:2));
                  fi_g = fi_g + D(i)*Beta_c1(i)*phi_alpha(sigma_norm(q_beta - Agents{1,i}(j,1:2)),d_prime(i),r_prime(i))*gradient_sigma_norm(q_beta - Agents{1,i}(j,1:2));
                  fi_d = fi_d + Beta_c2(i)*adj_matrix(q_beta - Agents{1,i}(j,1:2),d_prime(i))*(p_beta-Agents{1,i}(j,3:4));
              end
          end
          u_beta{1,i}(j,:) = fi_g + fi_d;
      end
  end
  
  u_beta = cell2mat(u_beta');
  
  u = u + u_beta;
end

agentForceVecX = u(:,1);
agentForceVecY = u(:,2);

%Repeat Mass (NAgents) number of times
mass_rep = repeat(NAgents,mass)';

out = [odeVec(2*TotAgent+1:4*TotAgent);agentForceVecX./mass_rep;agentForceVecY./mass_rep];

end

function out = repeat(r,x)
    t = r > 0;
    a = cumsum(r(t));
    b = zeros(1,a(end));
    b(a - r(t) + 1) = 1;
    x1 = x(t);
    out = x1(cumsum(b));
end

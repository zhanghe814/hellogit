function [out, components_size,components] = extract_connected_component(regulation_matrix)
% [regulation_matrix, kept_ind, components_size,components] = extract_connected_component(regulation_matrix)

regulation_matrix = double(regulation_matrix);
adj_matrix = (abs(regulation_matrix)+abs(regulation_matrix'))>0;
adj_matrix = adj_matrix - diag(diag(adj_matrix))+sparse(eye(size(adj_matrix)));

components=[];
is_assigned = zeros(size(regulation_matrix,1),1);
% is_assigned = full(sum(adj_matrix~=0)<=1); %%%% if we have this line, singleton nodes will not be considered as a component
while sum(is_assigned==0)~=0
%     sum(is_assigned==0);
%     [sum(is_assigned),length(is_assigned)]
    e = zeros(size(regulation_matrix,1),1);
    e(find(is_assigned==0,1))=1;
    while sum(e~=(adj_matrix*e>0))~=0
        e = (adj_matrix*e)>0;
    end
    components = [components,e];
    is_assigned(find(e==1))=1;
end

% components_size = sum(components);
% component_ind = find(components_size == max(components_size),1);
% 
% kept_ind = find(components(:,component_ind)==1);
% 
% regulation_matrix = regulation_matrix(kept_ind,kept_ind);
components_size = sum(components);
n=size(components_size,2);
for i=1:n
    kept_ind = find(components(:,i)==1);
    out(i).kept=kept_ind;
    out_matrix = regulation_matrix(kept_ind,kept_ind);
    out(i).matrix=out_matrix;
end

return
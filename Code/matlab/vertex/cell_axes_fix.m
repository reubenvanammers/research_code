function [length1, direction, chosen_vertices] = cell_axes_fix(l,C,V,fixed_vertices)

% vertex_list = C{l};
% nk = length(vertex_list);
% 
% A = -inf*ones(nk);
% 
% for i = 1:nk
%     for j = 1:i-1
%         A(i,j) = norm(V(vertex_list(i),:)-V(vertex_list(j),:));
%     end

% end
% 
% length1 = max(max(A));
% 
% [row, col] = find(A==length1,1); % pick one closest to x axis?
% direction = V(vertex_list(row),:)-V(vertex_list(col),:);
% chosen_vertices = [vertex_list(row) vertex_list(col)];
chosen_vertices = fixed_vertices;
direction = V(fixed_vertices(1),:)-V(fixed_vertices(2),:);
length1 = norm(direction);

if direction(1) < 0
    direction = -direction;
end

direction = direction/length1;
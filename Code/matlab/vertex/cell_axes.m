function [length1, direction, chosen_vertices] = cell_axes(l,C,V)
%Calculates axial direction based upon distance between furthest vertices.
vertex_list = C{l};
nk = length(vertex_list);

A = -inf*ones(nk);



for i = 1:nk
    for j = 1:i-1
        A(i,j) = norm(V(vertex_list(i),:)-V(vertex_list(j),:));
    end
end

length1 = max(max(A));

[row, col] = find(A==length1,1); % pick one closest to x axis?
direction = V(vertex_list(row),:)-V(vertex_list(col),:);
chosen_vertices = [vertex_list(row) vertex_list(col)];
if direction(1) < 0
    direction = -direction;
    chosen_vertices = fliplr(chosen_vertices);
end

direction = direction/length1;
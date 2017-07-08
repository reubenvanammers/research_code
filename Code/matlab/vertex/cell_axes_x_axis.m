function [length1, direction, chosen_vertices] = cell_axes_x_axis(l,C,V)
x_axis_priority = 0.1; %How much angles close to x axis will be preferred

vertex_list = C{l};
nk = length(vertex_list);

A = -inf*ones(nk);



for i = 1:nk
    for j = 1:i-1
        A(i,j) = norm(V(vertex_list(i),:)-V(vertex_list(j),:));
    end
end

length1 = max(max(A));

[rowvec, colvec] = find(A-length1<= x_axis_priority); % pick one closest to x axis?
possible_vectors = V(vertex_list(rowvec),:)-V(vertex_list(colvec),:);
vector_angle = abs(possible_vectors(:,2)./possible_vectors(:,1));
[~,index] = min(vector_angle);
row = rowvec(index);
col = colvec(index);
length1 = A(row,col);direction = V(vertex_list(row),:)-V(vertex_list(col),:);
chosen_vertices = [vertex_list(row) vertex_list(col)];

if direction(1) < 0
    direction = -direction;
end

direction = direction/length1;
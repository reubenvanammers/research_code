function [length1, direction, chosen_vertices] = cell_axes(l,C,V,stickyvertices)


stickiness_value = 0.1; %How resilient a matrix is to retaining old axes
%rather than choosing new ones; higher values mean will retain old values
%more. Must be greater than or equal to 0. 

x_axis_priority = 0.1; %How much angles close to x axis will be preferred

if nargin <4 
    stickyvertices = [];
end
vertex_list = C{l};
nk = length(vertex_list);

A = -inf*ones(nk);

for i = 1:nk
    for j = 1:i-1
        A(i,j) = norm(V(vertex_list(i),:)-V(vertex_list(j),:));
    end
end

length1 = max(max(A));
if stickyvertices ~= []
    if abs(A(stickyvertices(1),stickyvertices(2)) - length1) < stickiness_value;
        row = stickyvertices(1);
        col = stickyvertices(2);
        length1 = norm(V(stickyvertices(1),:),V(stickyvertices(2),:));
    end

else
    [rowvec, colvec] = find(A-length1<= x_axis_priority); % pick one closest to x axis?
    possible_vectors = V(vertex_list(rowvec),:)-V(vertex_list(colvec),:);
    vector_angle = possible_vectors(:,1)./possible_vectors(:,2);
    [~,index] = min(vector_angle);
    row = rowvec(index);
    col = colvec(index);
    length1 = A(row,col);
end
direction = V(vertex_list(row),:)-V(vertex_list(col),:);
chosen_vertices = [vertex_list(row) vertex_list(col)];
if direction(1) < 0
    direction = -direction;
end

direction = direction/length1;
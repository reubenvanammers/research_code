function d = d_lj(l,j,C,V); 
vertex_list = C{l};
%computes distance between j and next anticlockwise element of the cell
%represented by vertex_list, with local ordering j (represented by a list of logicals.)
%ie j = [0 0 0 1 0 0] represents the fourth vertex.
d = norm(V(vertex_list(logical(j)),:)-V(vertex_list(circshift(logical(j),1,2)),:));
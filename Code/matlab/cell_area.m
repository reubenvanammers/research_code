function area = cell_area(l,C,V);
vertex_list = C{l};
A = 0;
nk = length(vertex_list);
for j = 1:nk
    A = A + V(vertex_list(j),1)*V(vertex_list(mod(j,nk)+1),2)-V(vertex_list(mod(j,nk)+1),1)*V(vertex_list(j),2);
end
area = 0.5*abs(A);
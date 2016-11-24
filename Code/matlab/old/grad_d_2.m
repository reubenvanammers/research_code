function d_vec =grad_d_2(i,l,index,C,V)
%implements gradient of cell distance for use in vertex forces based on
%circumference. slower than grad_d_3
vertex_list = C{l};
%vertex i, cell l, index (ie, -1, 1), list of vertices locations V and  vertex_list Cell
%structures C{l}.
j=circshift((vertex_list==i),index,2);
i2 = circshift(j,1,2);
dlj = d_lj(l,j,C,V);
d_vec = 1./dlj.*(V(vertex_list(i2),:)-V(vertex_list(j),:));
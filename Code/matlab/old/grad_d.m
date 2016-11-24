function d_vec =grad_d(i,l,index,C,V)
%doesnt work use grad_d_2 or grad_d_3
vertex_list = C{l};
%vertex i, cell l, index (ie, -1, 1), list of vertices locations V and  vertex_list Cell
%structures C{l}.
j=circshift((vertex_list==i),index,2);

dlj = d_lj(l,j,C,V);
d_vec = 1./dlj.*(V(i,:)-V(vertex_list(j),:));
function d_vec = grad_d_tot(i,l,C,V)

len = length(C{l});
internal_vertex_number = find(i==C{l})-1;
im1 = C{l}(mod(internal_vertex_number -1,len)+1);
ip1 = C{l}(mod(internal_vertex_number + 1,len)+1);

V1 = V(i,:)-V(im1,:);
V1norm = V1./norm(V1);
V2 = V(i,:)-V(ip1,:);
V2norm = V2./norm(V2);
d_vec = V1norm +V2norm;
function d_vec = grad_d_3(i,l,index,C,V)
%implements gradient of cell distance for use in vertex forces based on
%circumference.

len = length(C{l});
internal_vertex_number = find(i==C{l})-1;
j2 = C{l}(mod(internal_vertex_number + index,len)+1);
i2 = C{l}(mod(internal_vertex_number + index + 1,len)+1);


% if index == -1;
%     j2 = neighbouring_cells{i}{l}(1);
%     i2 = neighbouring_cells{i}{l}(2);
% else
%     j2 = neighbouring_cells{i}{l}(2);
%     i2 = neighbouring_cells{i}{l}(3);
% end
dlj = norm(V(i2,:)-V(j2,:));
d_vec = 1./dlj.*(V(i2,:)-V(j2,:));
    
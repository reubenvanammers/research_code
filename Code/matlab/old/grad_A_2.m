function A_vec = grad_A_2(i,l,C,V)
%implements gradient of area for area based force in vertex model
% j = C{l}==i;
% jp1 = circshift(j,1,2);
% jm1 = circshift(j,-1,2);
len = length(C{l});
internal_vertex_number = find(i==C{l})-1;
vp1 = V(C{l}(mod(internal_vertex_number + 1,len)+1),:);
vm1 = V(C{l}(mod(internal_vertex_number - 1,len)+1),:);

% vp1 = V(C{l}(jp1),:);
% vm1 = V(C{l}(jm1),:);
A_vec = 0.5*[vp1(2)-vm1(2);vm1(1)-vp1(1)]';
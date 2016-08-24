function A_vec = grad_A_2(i,l,C,V)
j = C{l}==i;
jp1 = circshift(j,1,2);
jm1 = circshift(j,-1,2);
vp1 = V(C{l}(jp1),:);
vm1 = V(C{l}(jm1),:);
A_vec = 0.5*[vp1(2)-vm1(2);vm1(1)-vp1(1)]';
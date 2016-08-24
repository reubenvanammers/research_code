function A_vec = grad_A(i,l,C,V)
 A_vec = 0.5*[V(C{l}(circshift(C{l}==i,1,2)),:)-V(C{l}((circshift(C{l}==i,-1,2))))]*[0 -1; 1 0];

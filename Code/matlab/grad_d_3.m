function d_vec = grad_d_3(i,l,index,C,V,neighbouring_cells)
%implements gradient of cell distance for use in vertex forces based on
%circumference. 
if index == -1;
    j2 = neighbouring_cells{i}{l}(1);
    i2 = neighbouring_cells{i}{l}(2);
else
    j2 = neighbouring_cells{i}{l}(2);
    i2 = neighbouring_cells{i}{l}(3);
end
dlj = norm(V(i2,:)-V(j2,:));
d_vec = 1./dlj.*(V(i2,:)-V(j2,:));
    
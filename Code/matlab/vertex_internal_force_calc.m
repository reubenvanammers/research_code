function vertex_force = vertex_internal_force_calc(C,V,included_cell,lambda,beta,gamma,A_target,C_target)
%calculates the internal forces felt by vertices on connectivity and target
%parameters. Used for both real and reference cells. 
%A_target(i),C_target(i) should be vectors for each cell
N = length(V);
M = length(C);
vertex_force = zeros(N,2);
for  i= 1:N;%i is vertex number
    for l = included_cell{i};%l is cell number
        gradilm1 = grad_d_3(i,l,-1,C,V);
        gradil0 = grad_d_3(i,l,0,C,V);
        vertex_force(i,:) = vertex_force(i,:) +...
            2*lambda*(cell_area(l,C,V)-A_target(l))*grad_A_2(i,l,C,V)+...
            2*beta*(cell_circumference(l,C,V)-C_target(l))*(gradilm1-gradil0)+...
            gamma*gradilm1-gamma*gradil0;%due to definition subtract grad_d_3 instead of add
    end
end
vertex_force = -vertex_force;
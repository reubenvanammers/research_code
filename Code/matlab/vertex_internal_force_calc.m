function vertex_force = vertex_internal_force_calc(C,V,included_cell,lambda,beta,gamma,A_target,A_current,C_target,C_current)
%calculates the internal forces felt by vertices on connectivity and target
%parameters. Used for both real and reference cells. 
%A_target(i),C_target(i) should be vectors for each cell
N = length(V);
M = length(C);
vertex_force = zeros(N,2);
for  i= 1:N;%i is vertex number
    for l = included_cell{i};%l is cell number
        grad_d_vec = grad_d_tot(i,l,C,V);
        vertex_force(i,:) = vertex_force(i,:) +...
            2*lambda*(A_current(l)-A_target(l))*grad_A_2(i,l,C,V)+...
            2*beta*(C_current(l)-C_target(l))*(grad_d_vec)+...
            gamma*(grad_d_vec);
    end
end
vertex_force = -vertex_force;
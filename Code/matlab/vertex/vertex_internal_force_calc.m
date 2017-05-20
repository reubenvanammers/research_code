function vertex_force = vertex_internal_force_calc(C,V,included_cell,lambda,beta,gamma,A_target,A_current,C_target,C_current)
%calculates the internal forces felt by vertices on connectivity and target
%parameters. Used for both real and reference cells. 
%A_target(i),C_target(i) should be vectors for each cell
N = length(V);
M = length(C);
vertex_force = zeros(N,2);
for  i= 1:N;%i is vertex number
    for l = included_cell{i};%l is cell number
        len = length(C{l});
        internal_vertex_number = find(i==C{l})-1;
        v = V(i,:);
        vp1 = V(C{l}(mod(internal_vertex_number + 1,len)+1),:);
        vm1 = V(C{l}(mod(internal_vertex_number - 1,len)+1),:);
        
        
        grad_A_vec = 0.5*[vp1(2)-vm1(2);vm1(1)-vp1(1)]';
        
        V1 = v-vm1;
        V2 = v-vp1;
        V1norm = V1./norm(V1);
        V2norm = V2./norm(V2);
        
        grad_d_vec = V1norm +V2norm;

        vertex_force(i,:) = vertex_force(i,:) +...
            2*lambda*(A_current(l)-A_target(l))*grad_A_vec+...
            2*beta*(C_current(l)-C_target(l))*(grad_d_vec)+...
            gamma*(grad_d_vec);
    end
end
vertex_force = -vertex_force;
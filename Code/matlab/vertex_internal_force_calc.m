function vertex_force = vertex_internal_force_calc(connectivitylist,C,V,lambda,beta,gamma,A_target,C_target)
%A_target(i),C_target(i) should be vectors for each cell
N = length(V);
M = length(C);
vertex_force = zeros(N,2);
for  i= 2:N;
    for l = 1:M;
        if connectivitylist(l,i) ==1;
            vertex_force(i,:) = vertex_force(i,:) +...
                2*lambda*(cell_area(l,C,V)-A_target(l))*grad_A_2(i,l,C,V)+...
                2*beta*(cell_circumference(l,C,V)-C_target(l))*(grad_d_2(i,l,-1,C,V)-grad_d_2(i,l,0,C,V))+...
                gamma*grad_d_2(i,l,-1,C,V)-gamma*grad_d_2(i,l,0,C,V);%due to definition subtract grad_d_2 instead of add
        end
    end
end
vertex_force = -vertex_force;
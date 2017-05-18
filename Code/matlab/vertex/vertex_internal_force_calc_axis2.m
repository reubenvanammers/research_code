function vertex_force = vertex_internal_force_calc_axis2(C,V,included_cell,lambda,beta,gamma,A_target,A_current,C_target,C_current,axis_target,axis_current)
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
            0*gamma*(grad_d_vec); %removing this from the equation 
    end
end

for l = 1:M
        l_target = axis_target{l}{1};
        direction_target = axis_target{l}{2};
        
        
        l_current = axis_current{l}{1};
        v1 = V(axis_current{l}{3}(1),:);
        v2 = V(axis_current{l}{3}(2),:);
        
        vlen = norm(v1-v2);
        
        length_force = 2*(l_current-l_target)*(v1-v2)/vlen;
        
        angle_force = ((v1(1)-v2(1))*direction_target(1)+(y1-y2)*direction_target(2))*(v1-v2)/(-vlen^3) +...
            direction_target;
        
        total_force = length_force+angle_force;
        
        vertex_force(axis_current{l}{3}(1),:) = vertex_force(axis_current{l}{3}(1),:) + gamma*total_force;
        vertex_force(axis_current{l}{3}(2),:) = vertex_force(axis_current{l}{3}(2),:) - gamma*total_force;
end
       
        
        
vertex_force = -vertex_force;
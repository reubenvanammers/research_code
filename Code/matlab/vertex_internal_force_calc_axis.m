function vertex_force = vertex_internal_force_calc_axis(C,V,included_cell,lambda,beta,gamma,A_target,A_current,C_target,C_current,axis_target,Ixx,Iyy,Ixy,Cx,Cy)
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
        
        z = [v(1)*vp1(2)-vp1(1)*v(2);vm1(1)*v(2)-v(1)*vm1(2)]; %Commonly used value
        
        grad_A_vec = 0.5*[vp1(2)-vm1(2);vm1(1)-vp1(1)]';
        
        V1 = v-vm1;
        V2 = v-vp1;
        V1norm = V1./norm(V1);
        V2norm = V2./norm(V2);
        
        grad_Ixx = 1/12*[
            (v(2)^2+v(2)*vp1(2)+vp1(2)^2)*vp1(2)+(vm1(2)^2+vm1(2)*v(2)+v(2)^2)*(-vm1(2));...
            (2*v(2)+vp1(2))*z(1)+(v(2)^2+v(2)*vp1(2)+vp1(2)^2)*(-vp1(1)) + ...
            (2*v(2)+vm1(2))*z(2)+(vm1(2)^2+v(2)*vm1(2)+v(2)^2)*(vm1(1))];
        
        grad_Iyy = 1/12*[
            (2*v(1)+vp1(1))*z(1)+(v(1)^2+v(1)*vp1(1)+vp1(1)^2)*(vp1(2))+ ...
            (2*v(1)+vm1(1))*z(2)+(v(1)^2+v(1)*vm1(1)+vm1(1)^2)*(-vm1(2));
            (v(1)^2+v(1)*vp1(1)+vp1(1)^2)*(-vp1(1))+(vm1(1)^2+vm1(1)*v(1)+v(1)^2)*vm1(1)];
        
        grad_Ixy = 1/24*[
            (2*v(2)+vp1(2))*z(1)+(v(1)*vp1(2)+2*v(1)*v(2)+2*vp1(1)*vp1(2)+vp1(1)*v(2))*vp1(2)+ ...
            (2*v(2)+vm1(2))*z(2)+(vm1(1)*v(2)+2*vm1(1)*vm1(2)+2*v(1)*v(2)+v(1)*vm1(2))*(-vm1(2));
            (2*v(1)+vp1(1))*z(1)+(v(1)*vp1(2)+2*v(1)*v(2)+2*vp1(1)*vp1(2)+vp1(1)*v(2))*(-vp1(1)) + ...
            (2*v(1)+vm1(1))*z(2)+(vm1(1)*v(2)+2*vm1(1)*vm1(2)+2*v(1)*v(2)+v(1)*vm1(2))*(vm1(1))];
        
        
        grad_P = grad_Iyy-grad_Ixx+(Ixx(l)*grad_ixx+Iyy(l)*grad_Iyy-4*Ixy(l)*grad_Ixy-Ixx(l)*grad_Iyy-Iyy(l)*grad_Ixx)/...
            ((Ixx(l)^2+Iyy(l)^2-4*Ixy(l)^2-2*Ixx(l)*Iyy(l))^0.5);
        
        P = Iyy(l)-Ixx(l)+(Ixx(l)^2+Iyy(l)^2-4*Ixy(l)^2-2*Ixx(l)*Iyy(l))^0.5;
        
        grad_d_vec = V1norm +V2norm;
        
        axis_force = (P*grad_P+4*Ixy(l)*grad_Ixy)*(P*axis_target(1)+2*Ixy(l)*axis_target(2))/((P^2+4*Ixy)^1.5)-...
            (((P^2+4*Ixy)^(-0.5))*(grad_P*axis_target(1)+grad_Ixy*axis_target(2)));

        vertex_force(i,:) = vertex_force(i,:) +...
            2*lambda*(A_current(l)-A_target(l))*grad_A_vec+...
            2*beta*(C_current(l)-C_target(l))*(grad_d_vec)+...
            gamma*(grad_d_vec);
    end
end
vertex_force = -vertex_force;
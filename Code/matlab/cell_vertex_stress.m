function dxdt = cell_vertex_stress(t,x);

global C connectivitylist F N A0 C0 lambda beta gamma M
dxdt = zeros(2*N,1);
vertex_force = zeros(N,2);
[V,~] = matricize([x;x]);
for  i= 2:N;
    for l = 1:M;
        if connectivitylist(l,i) ==1;
            vertex_force(i,:) = vertex_force(i,:) +...
                2*lambda*(cell_area(l,C,V)-A0)*grad_A_2(i,l,C,V)+...
                2*beta*(cell_circumference(l,C,V)-C0)*(grad_d_2(i,l,-1,C,V)-grad_d_2(i,l,0,C,V))+...
                gamma*grad_d_2(i,l,-1,C,V)-gamma*grad_d_2(i,l,0,C,V);%due to definition subtract grad_d_2 instead of add
        end
    end
end

dxdt = columnize(-vertex_force,-vertex_force);
dxdt = dxdt(1:2*N);
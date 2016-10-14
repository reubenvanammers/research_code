function dxdt = cell_vertex_stress(t,x);
%used in simulation without remodelling.
global C F N A0_vec C0_vec lambda beta gamma M included_cell
dxdt = zeros(2*N,1);
vertex_force = zeros(N,2);
[V,~] = matricize([x;x]);
t
real_cell_areas = zeros(1,M);
real_cell_circumferences = zeros(1,M);

for i = 1:M
    real_cell_areas(i) = cell_area(i,C,V);
    real_cell_circumferences(i) = cell_circumference(i,C,V);
end

force = vertex_internal_force_calc(C,V,included_cell,lambda,beta,gamma,A0_vec,real_cell_areas,C0_vec,real_cell_circumferences);



dxdt = columnize(force,force);
dxdt = dxdt(1:2*N);
function dxdt = cell_vertex_stress_reference(t,x);

global C connectivitylist F N A0_vec C0_vec lambda beta gamma M alpha

[V,V_ref] = matricize([x]);
real_cell_areas = zeros(1,M);
real_cell_circumferences = zeros(1,M);
reference_cell_areas = zeros(1,M);
reference_cell_circumferences = zeros(1,M);
for i = 1:M
    real_cell_areas(i) = cell_area(i,C,V);
    real_cell_circumferences(i) = cell_circumference(i,C,V);
    reference_cell_areas(i) = cell_area(i,C,V_ref);
    reference_cell_circumferences(i) = cell_circumference(i,C,V_ref);
end

real_force = vertex_internal_force_calc(connectivitylist,C,V,lambda,beta,gamma,reference_cell_areas,reference_cell_circumferences);
follow_force = (1-alpha).*vertex_internal_force_calc(connectivitylist,C,V_ref,lambda,beta,gamma,real_cell_areas,real_cell_circumferences);
fix_force = alpha*vertex_internal_force_calc(connectivitylist,C,V_ref,lambda,beta,gamma,A0_vec,C0_vec);
t
dxdt = columnize(real_force,follow_force+fix_force);
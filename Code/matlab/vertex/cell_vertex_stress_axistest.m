function dxdt = cell_vertex_stress_axistest(t,x)
%used in simulation without remodelling.
global C N A0_vec C0_vec lambda beta gamma delta M included_cell V axis_target fixed_vertices
t
dxdt = zeros(2*N,1);
vertex_force = zeros(N,2);
[V,~] = matricize([x;x]);

real_cell_areas = zeros(1,M);
real_cell_circumferences = zeros(1,M);
axis_current = cell(1,M);

for i = 1:M
    real_cell_areas(i) = cell_area(i,C,V);
    real_cell_circumferences(i) = cell_circumference(i,C,V);
    [len,direction,vertices] = cell_axes_fix(i,C,V,fixed_vertices{i}); % give stickingess?
    axis_current{i} = {len,direction,vertices};
end

force = vertex_internal_force_calc_axis2(C,V,included_cell,lambda,beta,gamma,delta,A0_vec,real_cell_areas,C0_vec,real_cell_circumferences,axis_target,axis_current);



dxdt = columnize(force,force);
dxdt = dxdt(1:2*N);
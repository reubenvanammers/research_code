global V C connectivitylist F N A0_vec C0_vec lambda beta gamma M included_cell 
%vertex model without remodelling
sidelength = 1/sqrt(3)*1.5;
A0=sqrt(27)/2*(sidelength.^2);
C0 = 6*sidelength;
lambda = 1;
beta = 1;
gamma = 0;
[V,C,connectivitylist] = hexgrid_voronoi(5,5);
V_init = V;
N= length(V);
M = length(C);
A0_vec = ones(1,M)*A0;
C0_vec = ones(1,M)*C0;
V(3,1) = V(3,1)+0;
V_ref = V;
V_vec = columnize(V,V_ref);
V_vec = V_vec(2*N+1:end);%ignore reference cells for now
F=0;
tend = 1000;
included_cell = cell_inclusion(V,C);

options = odeset('RelTol',1e-3,'AbsTol',1e-6);
[Time,Y] = ode15s(@cell_vertex_stress,0:0.2:tend,V_vec,options);
final_hex = Y(end,:)'
[V,~] = matricize([final_hex;final_hex])

cell_areas = zeros(1,M);
cell_circumferences = zeros(1,M);
for i = 1:M
    cell_areas(i) = cell_area(i,C,V);
    cell_circumferences(i) = cell_circumference(i,C,V);
end
area_diff = norm(cell_areas-A0)
circumference_diff = norm(cell_circumferences-C0)
hex_vis(Time,Y,C)

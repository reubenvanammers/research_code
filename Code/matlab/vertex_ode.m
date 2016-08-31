global V C connectivitylist F N A0 C0 lambda beta gamma M

A0=sqrt(3)/2;
C0 = 2*sqrt(pi*A0);
lambda = 0;
beta = 1;
gamma = 0;
[V,C,connectivitylist] = hexgrid_voronoi();
V_init = V;
N= length(V);
M = length(C);
V(3,1) = V(3,1)+0;
ref_V = V;
V_vec = columnize(V,ref_V);
V_vec = V_vec(1:2*N);%ignore reference cells for now
F=0;
tend = 100;


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

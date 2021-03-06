function [Time,Y,C2] = vertex_ode_reference_nocirc(lambda0,beta0,gamma0,alpha0,eta0,T0,tend)
%implements vertex model with remodelling
%tries to match real(reference) circumference with reference (real) area
%instead of circumference, or vice versa
global C N A0_vec C0_vec lambda beta gamma M alpha circ_area_conversion external_force
global t_rec C_rec A_rec T fixlist movelist eta restoring_rec counter included_cell
sidelength = 1/sqrt(3);
A0=sqrt(27)/2*(sidelength.^2);

C0 = 6*sidelength;

circ_area_conversion = 2;
%this gives the type of interplay between the circumference and the area
%see cell_vertex_stress_reference_nocirc for details

lambda = lambda0;beta=beta0;gamma=gamma0;alpha=alpha0;T=T0;
eta = eta0;
[V,C] = hexgrid_voronoi([10 10]);
included_cell = cell_inclusion(V,C);
external_force = 1;
N= length(V);
M = length(C);
A0_vec = ones(1,M)*A0;
C0_vec = ones(1,M)*C0;
V(3,1) = V(3,1)+0;
ref_V = V;
V_vec = columnize(V,ref_V);


m = min(V(:,1));
fixlist = V(:,1) <m+0.4;
m = max(V(:,1));
movelist =  V(:,1) >m-0.4;%plus minus 0.1 is for minor discrepancies


t_rec = linspace(-T,0)';%sets up averaging vector for each edge
counter = 1;
C_rec = C0*ones(100,M);
A_rec = A0*ones(100,M);
restoring_rec = [];

options = odeset('RelTol',1e-5,'AbsTol',1e-8,'Events',@intersect_event);
[Time,Y] = ode15s(@cell_vertex_stress_reference_nocirc,0:0.001:tend,V_vec,options);
%final_hex = Y(end,:)';
C2 = C;
%hex_vis_2(Time,Y,C);


% cell_areas = zeros(1,M);
% cell_circumferences = zeros(1,M);
% for i = 1:M
%     cell_areas(i) = cell_area(i,C,V);
%     cell_circumferences(i) = cell_circumference(i,C,V);
% end
% area_diff = norm(cell_areas-A0);
% circumference_diff = norm(cell_circumferences-C0);

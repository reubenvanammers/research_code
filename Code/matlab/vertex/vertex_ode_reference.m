function [Time,Y,C2] = vertex_ode_reference(lambda0,beta0,gamma0,alpha0,eta0,T0,tend)
%implements vertex model with remodelling
global C N A0_vec C0_vec lambda beta gamma M alpha external_force
global t_rec C_rec A_rec T fixlist movelist eta restoring_rec counter included_cell
sidelength = 1/sqrt(3);
A0=sqrt(27)/2*(sidelength.^2);
%C0 = 2*sqrt(pi*A0);
C0 = 6*sidelength;

external_force = 7.57;

lambda = lambda0;beta=beta0;gamma=gamma0;alpha=alpha0;T=T0;
eta = eta0;
[V,C] = hexgrid_voronoi([10,10]);

included_cell = cell_inclusion(V,C);
N= length(V);
M = length(C);
A0_vec = ones(1,M)*A0;
C0_vec = ones(1,M)*C0;
V(3,1) = V(3,1)+0;
ref_V = V;
V_vec = columnize(V,ref_V);


m = min(V(:,1));
fixlist = V(:,1) <m+0.1;
m = max(V(:,1));
movelist =  V(:,1) >m-0.1;%plus minus 0.1 is for minor discrepancies


t_rec = linspace(-T,0)';%sets up averaging vector for each edge
counter = 1;
C_rec = C0*ones(100,M);
A_rec = A0*ones(100,M);
restoring_rec = [];

options = odeset('RelTol',1e-5,'AbsTol',1e-8);
[Time,Y] = ode15s(@cell_vertex_stress_reference,0:0.2:tend,V_vec,options);

C2 =C;

% final_hex = Y(end,:)';
% [V,V_ref] = matricize(final_hex);
%hex_vis_2(Time,Y,C);


% cell_areas = zeros(1,M);
% cell_circumferences = zeros(1,M);
% for i = 1:M
%     cell_areas(i) = cell_area(i,C,V);
%     cell_circumferences(i) = cell_circumference(i,C,V);
% end
% area_diff = norm(cell_areas-A0);
% circumference_diff = norm(cell_circumferences-C0);

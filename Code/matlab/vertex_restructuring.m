function [Time,Y,cell_history2,cell_t_history2] = vertex_restructuring(lambda0,beta0,gamma0,alpha0,eta0,T0,tend)
%implements vertex model with remodelling
global C F N A0_vec C0_vec lambda beta gamma M alpha cell_history cell_t_history
global t_rec C_rec A_rec T fixlist movelist eta restoring_rec counter included_cell external_force
sidelength = 1/sqrt(3);
A0=sqrt(27)/2*(sidelength.^2);
C0 = 6*sidelength;

lambda = lambda0;beta=beta0;gamma=gamma0;alpha=alpha0;T=T0;
eta = eta0;
[V,C] = hexgrid_voronoi();
external_force = 0.2;
included_cell = cell_inclusion(V,C);
N= length(V);
M = length(C);
A0_vec = ones(1,M)*A0;
C0_vec = ones(1,M)*C0;z
V(3,1) = V(3,1)+0;
ref_V = V;
V_vec = columnize(V,ref_V);


F=@stress_force_sync;
m = min(V(:,1));
fixlist = V(:,1) <m+0.1;
m = max(V(:,1));
movelist =  V(:,1) >m-0.1;%plus minus 0.1 is for minor discrepancies


t_rec = linspace(-T,0)';%sets up averaging vector for each edge
counter = 1;
C_rec = C0*ones(100,M);
A_rec = A0*ones(100,M);
restoring_rec = [];

options = odeset('RelTol',1e-3,'AbsTol',1e-6);
[Time,Y] = cell_intermediate_restructuring(@cell_vertex_stress_reference,tend,V_vec,options);
final_hex = Y(end,:)';
cell_history2 =cell_history;
cell_t_history2 = cell_t_history;
%hex_vis_2(Time,Y,C);



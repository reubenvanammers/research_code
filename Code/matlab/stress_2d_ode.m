function [Time,Y]=stress_2d_ode(alpha0,gamma0,T0,tend)
%Implements remodelling in a cell centre cell centre spring based model.
global gamma alpha s0 F N M E r_rec t_rec T fixlist movelist vertex_matrix_1 vertex_matrix_2 edge_matrix restoring_rec restoring_t_rec
gamma = gamma0;alpha =alpha0;T = T0;
%Model Parameters
%gamma = 1;%speed of reference cell remodelling
%T= 20; %time over which length is remembered
%alpha = 0.5; % proportion of reference cell remodelling based on default length
s0 = 1; %un-stretch length of cells 
restoring_rec = [];
restoring_t_rec = [];
F=@stress_force_sync;



[P,E,Tri] = tri_2D_Hex2;
m = min(P(:,1));
fixlist = P(:,1) ==m;
m = max(P(:,1));
movelist =  P(:,1) ==m;


N= length(P);
M = length(E);

[vertex_matrix_1,vertex_matrix_2,edge_matrix] = edge_matrix_create(E,N);


r_rec = s0*ones(100,M);
t_rec = linspace(-T,0)';%sets up averaging vector for each edge
ref_P = P;
tot_P = columnize(P,ref_P);
%essentially have 4 lists of data stacked in one column vector:
% real x values, real y values, reference x values, reference y values
options = odeset('RelTol',1e-5,'AbsTol',1e-8);
[Time,Y] = ode15s(@cell_forces_stress_vector,0:0.2:tend,tot_P,options);

%tri_vis(Time,Y,Tri)%visualizes system

end



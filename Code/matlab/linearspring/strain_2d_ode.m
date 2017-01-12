function [Time,Y,Tri2,stress_rec2]=strain_2d_ode(alpha0,eta0,T0,tend,strainfunc,t_strain_end2,gridsize)
%Implements remodelling in a cell centre cell centre spring based model.
%Needs external strain function value to match that strain. Optional
%argument when to stop the external strain.
global eta alpha s0 N M E r_rec t_rec T fixlist movelist vertex_matrix_1 k initial_length t_strain_end
global vertex_matrix_2 edge_matrix restoring_rec restoring_t_rec counter strain_function stress_rec
eta = eta0;alpha =alpha0;T = T0;
%Model Parameters
%eta = 1;%speed of reference cell remodelling
%T= 20; %time over which length is remembered
%alpha = 0.5; % proportion of reference cell remodelling based on default length
s0 = 1; %un-stretch length of cells 
k=1; %spring constant
restoring_rec = [];
stress_rec = [];
restoring_t_rec = [];

if nargin < 7
    gridsize = [7,8];
end
if nargin < 6
    t_strain_end = Inf;
else
    t_strain_end = t_strain_end2;
end

strain_function = strainfunc;
[P,E,Tri] = tri_2D_Hex2(gridsize);
initial_min = min(P(:,1));
fixlist = P(:,1) ==initial_min;
initial_max = max(P(:,1));
movelist =  P(:,1) ==initial_max;
initial_length = initial_max-initial_min;

N= length(P);
M = length(E);

[vertex_matrix_1,vertex_matrix_2,edge_matrix] = edge_matrix_create(E,N);


r_rec = s0*ones(100,M);
t_rec = linspace(-T,0)';%sets up averaging vector for each edge
counter = 1;
ref_P = P;
tot_P = columnize(P,ref_P);
%essentially have 4 lists of data stacked in one column vector:
% real x values, real y values, reference x values, reference y values
options = odeset('RelTol',1e-5,'AbsTol',1e-8);
[Time,Y] = ode15s(@cell_forces_strain_vector,0:0.2:tend,tot_P,options);
Tri2 = Tri;
%tri_vis(Time,Y,Tri)%visualizes system


l = length(t_rec)-length(stress_rec);
t_rec2 = t_rec(l+1:end);
[t_rec2,ia,~] = unique(t_rec2);
stress_rec2 = stress_rec(ia);%deletes duplicate time entries for interpolation
stress_rec2 = interp1(t_rec2,stress_rec2,Time);
end



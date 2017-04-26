function [Time,Y,Tri2,flag]=stress_2d_ode(alpha0,eta0,T0,tend,gridsize,ext_force,maxstrain)
%Implements remodelling in a cell centre cell centre spring based model.
%Has external stress force applied to rightmost edge of cells.
global eta alpha s0 N M E r_rec t_rec T fixlist movelist vertex_matrix_1 k 
global vertex_matrix_2 edge_matrix restoring_rec restoring_t_rec counter external_force
global maxlength strainflag
eta = eta0;alpha =alpha0;T = T0;
strainflag = false;

if nargin < 5
    gridsize = [7,8]; %default size of monolayer
end
if nargin < 6
    external_force = 0.2;
else
    external_force = ext_force;
end

%Model Parameters
%eta = 1;%speed of reference cell remodelling
%T= 20; %time over which length is remembered
%alpha = 0.5; % proportion of reference cell remodelling based on default length
s0 = 1; %un-stretch length of cells 
k=1; %spring constant
restoring_rec = [];
restoring_t_rec = [];


[P,E,Tri] = tri_2D_Hex2(gridsize);
initial_min = min(P(:,1));
fixlist = P(:,1) ==initial_min;
initial_max = max(P(:,1));
movelist =  P(:,1) ==initial_max;
if nargin ==7
    maxlength = maxstrain*(initial_max-initial_min);
end


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
options = odeset('RelTol',1e-5,'AbsTol',1e-8,'Events',@stress_event);
[Time,Y] = ode15s(@cell_forces_stress_vector_maxstrain,0:0.2:tend,tot_P,options);
Tri2 = Tri;
%tri_vis(Time,Y,Tri)%visualizes system
flag = strainflag;
end



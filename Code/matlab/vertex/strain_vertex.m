function [Time,Y,C2,stress_rec2,t_rec2,stress_index] = strain_vertex(lambda0,beta0,gamma0,alpha0,eta0,T0,tend,strainfunc,gridsize,t_ramp_end,t_strain_end2)
%implements vertex model with remodelling
%tries to match real(reference) circumference with reference (real) area
%instead of circumference. 
global C N A0_vec C0_vec lambda beta gamma M alpha circ_area_conversion  
global t_rec C_rec A_rec T fixlist movelist eta restoring_rec counter included_cell strain_function initial_length t_strain_end
global stress_rec restoring_t_rec
sidelength = 1/sqrt(3);
A0=sqrt(27)/2*(sidelength.^2);
circ_area_conversion = 1;



if nargin < 11
    t_strain_end = Inf;
else
    t_strain_end = t_strain_end2;
end
if nargin < 10
    t_ramp_end = tend;
end
if nargin < 9
    gridsize = [7,8]; %default size of monolayer
end


strain_function = strainfunc;

C0 = 6*sidelength;


%this gives the type of interplay between the circumference and the area
%see cell_vertex_stress_reference_nocirc for details

lambda = lambda0;beta=beta0;gamma=gamma0;alpha=alpha0;T=T0;
eta = eta0;
[V,C] = hexgrid_voronoi(gridsize);
included_cell = cell_inclusion(V,C);
%external_force = 0.2;
N= length(V);
M = length(C);
A0_vec = ones(1,M)*A0;
C0_vec = ones(1,M)*C0;
V(3,1) = V(3,1)+0;
ref_V = V;
V_vec = columnize(V,ref_V);


initial_min = min(V(:,1));
fixlist = V(:,1) <initial_min+0.4;
initial_max = max(V(:,1));
movelist =  V(:,1) >initial_max-0.4;%plus minus 0.4 is for first two layers of cells
initial_length = initial_max-initial_min;

t_rec = linspace(-T,0)';%sets up averaging vector for each edge
counter = 1;
C_rec = C0*ones(100,M);
A_rec = A0*ones(100,M);
restoring_rec = [];
stress_rec = [];
restoring_t_rec = [];
options = odeset('RelTol',1e-5,'AbsTol',1e-8);
[Time,Y] = ode15s(@cell_vertex_strain_reference_nocirc,0:0.2:t_ramp_end,V_vec,options);
stress_index = length(stress_rec);
V_vec = Y(end,:)';
if t_ramp_end ~= tend
    [Time2,Y2] = ode15s(@cell_vertex_strain_reference_nocirc,t_ramp_end:tend,V_vec,options);
    %final_hex = Y(end)
    Time = [Time; Time2(2:end)];
    Y = [Y; Y2(2:end,:)];
end;
C2 = C;

l = length(t_rec)-length(stress_rec);
length(t_rec);
length(stress_rec);
t_rec2 = t_rec(l+1:end);
stress_rec2 = stress_rec;
%hex_vis_2(Time,Y,C);


% cell_areas = zeros(1,M);
% cell_circumferences = zeros(1,M);
% for i = 1:M
%     cell_areas(i) = cell_area(i,C,V);
%     cell_circumferences(i) = cell_circumference(i,C,V);
% end
% area_diff = norm(cell_areas-A0);
% circumference_diff = norm(cell_circumferences-C0);

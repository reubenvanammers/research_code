function [Time,Y,C2,stress_rec2,t_rec2,stress_index] = strain_vertex_axis(lambda0,beta0,gamma0,delta0,alpha0,eta0,T0,tend,gridsize,strainfunc,t_ramp_end,t_strain_end2)
%implements vertex model with remodelling
%tries to match real(reference) circumference with reference (real) area
%instead of circumference. 
global C N A0_vec C0_vec lambda beta gamma M alpha circ_area_conversion  
global t_rec C_rec A_rec T fixlist movelist eta restoring_rec counter included_cell strain_function initial_length t_strain_end
global stress_rec restoring_t_rec
global delta axis_0 fixed_vertices length_rec angle_x_proj_rec angle_y_proj_rec

sidelength = 1/sqrt(3);
A0=sqrt(27)/2*(sidelength.^2);
circ_area_conversion = 1;



if nargin < 12
    t_strain_end = Inf;
else
    t_strain_end = t_strain_end2;
end

if nargin < 9
    gridsize = [7,8]; %default size of monolayer
end


strain_function = strainfunc;

C0 = 6*sidelength;


%this gives the type of interplay between the circumference and the area
%see cell_vertex_stress_reference_nocirc for details

lambda = lambda0;beta=beta0;gamma=gamma0;alpha=alpha0;T=T0; delta = delta0;
eta = eta0;
[V,C] = hexgrid_voronoi(gridsize);
V(:,1) = V(:,1)*1.0000001;%slightly stretches x direction so that preferred direction is along x axis
included_cell = cell_inclusion(V,C);
%external_force = 0.2;
N= length(V);
M = length(C);
A0_vec = ones(1,M)*A0;
C0_vec = ones(1,M)*C0;
V(3,1) = V(3,1)+0;
ref_V = V;
V_vec = columnize(V,ref_V);


for i = 1:M
    %x = [1 0.5];
    %axis_target{i} = {1.9*sidelength,x/norm(x),[1 2]}; %length,direction, and some dummy axis values
    %(last entry is used to keep track of which vertices is used in length
    %and angle calculations, but is currently unused without a reference
    %state
    %[~,~,fixed_vertices{i}] = cell_axes(i,C,V); %fixing vertices to calculate from simplifies matters, 
    % and prevents solver from stalling when going below default target
    % length
    [base_cell_length,base_cell_direction,fixed_cell_vertices] = cell_axes(i,C,V);
    axis_0{i} = {base_cell_length base_cell_direction fixed_cell_vertices};
    fixed_vertices{i} = fixed_cell_vertices;
end



initial_min = min(V(:,1));
fixlist = V(:,1) <initial_min+0.1;
initial_max = max(V(:,1));
movelist =  V(:,1) >initial_max-0.1;%plus minus 0.1 is for minor discrepancies
initial_length = initial_max-initial_min;

t_rec = linspace(-T,0)';%sets up averaging vector for each edge
counter = 1;
C_rec = C0*ones(100,M);
A_rec = A0*ones(100,M);
restoring_rec = [];
stress_rec = [];
restoring_t_rec = [];

length_rec = ones(100,M);
angle_x_proj_rec = ones(100,M);
angle_y_proj_rec = ones(100,M);
for i = 1:M
    length_rec(:,i) = ones(100,1)*axis_0{i}{1};
    angle_x_proj_rec(:,i) = ones(100,1)*axis_0{i}{2}(1);
    angle_y_proj_rec(:,i) = ones(100,1)*axis_0{i}{2}(2);
end

options = odeset('RelTol',1e-5,'AbsTol',1e-8);
[Time,Y] = ode15s(@cell_vertex_strain_axis_reference,0:0.2:t_ramp_end,V_vec,options);
stress_index = length(stress_rec);
V_vec = Y(end,:)';
[Time2,Y2] = ode15s(@cell_vertex_strain_axis_reference,t_ramp_end:0.2:tend,V_vec,options);
%final_hex = Y(end
Time = [Time; Time2(2:end)];
Y = [Y; Y2(2:end,:)];

C2 = C;

l = length(t_rec)-length(stress_rec);
t_rec;
l;
length(t_rec)
length(stress_rec)
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

function [Time,Y]=ode_2D
close all
global gamma alpha s0 F N M E r_rec t_rec T
gamma = 1;
T= 10;
alpha = 0.2;
s0 = 1;
F=@ct_force;
tend = 50;


[P,E,Tri] = tri_2D_Hex;
N= length(P);
M = length(E);
r_rec = zeros(1,M);
t_rec = -T;%sets up averaging vector for each edge
ref_P = P;
tot_P = columnize(P,ref_P);
%essentially have 4 lists of data stacked in one column vector:
% real x values, real y values, reference x values, reference y values

[Time,Y] = ode15s(@cell_forces,[0:0.2:tend],tot_P);

tri_vis(Time,Y,Tri)%visualizes system

end



function [Time,Y]=ode_2D
%implements tension/compression test with remodelling.
close all
global gamma alpha s0 F N M E r_rec t_rec T
%Model Parameters
gamma = 2;%speed of reference cell remodelling
T= 10; %time over which length is remembered
alpha = 0.4; % proportion of reference cell remodelling based on default length
s0 = 1; %un-stretch length of cells 

F=@ct_force;
tend = 70;


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



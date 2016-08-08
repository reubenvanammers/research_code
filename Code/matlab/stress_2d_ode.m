function [Time,Y]=stress_2d_ode(alpha0,gamma0,T0)
%close all

global gamma alpha s0 F N M E r_rec t_rec T fixlist movelist
gamma = gamma0;alpha =alpha0;T = T0;
%Model Parameters
%gamma = 1;%speed of reference cell remodelling
%T= 20; %time over which length is remembered
%alpha = 0.5; % proportion of reference cell remodelling based on default length
s0 = 1; %un-stretch length of cells 

F=@stress_force;
tend = 1000;


[P,E,Tri] = tri_2D_Hex2;
m = min(P(:,1));
fixlist = P(:,1) ==m;
m = max(P(:,1));
movelist =  P(:,1) ==m;


N= length(P);
M = length(E);
r_rec = s0*ones(100,M);
t_rec = linspace(-T,0)';%sets up averaging vector for each edge
ref_P = P;
tot_P = columnize(P,ref_P);
%essentially have 4 lists of data stacked in one column vector:
% real x values, real y values, reference x values, reference y values

[Time,Y] = ode113(@cell_forces_stress,[0 tend],tot_P);

%tri_vis(Time,Y,Tri)%visualizes system

end



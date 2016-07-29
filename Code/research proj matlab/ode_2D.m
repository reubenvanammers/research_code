function [Time,Y]=ode_2D
close all
global gamma alpha s0 F N M E r_rec t_rec T
gamma = 3;
T= 1;
alpha = 1;
s0 = 0.75;
F=0;
tend = 20;


[P,E,Tri] = tri_2D;
N= length(P);
M = length(E);
r_rec = zeros(1,M);
t_rec = -T;%sets up averaging vector for each edge
P =P(:);
ref_P = P;
tot_P = [P;ref_P];
%essentially have 4 lists of data stacked in one column vector:
% real x values, real y values, reference x values, reference y values

[Time,Y] = ode15s(@cell_forces,[0:0.2:tend],tot_P);
for i = 1:length(Time);
    final_points = Y(i,:);
    real_points = reshape(final_points(1:2*N),[],2);
    TR = triangulation(Tri,real_points);
    triplot(TR);
    pause(0.1)
end

for i = 1:length(Time);
    final_points = Y(i,:);
    ref_points = reshape(final_points(2*N+1:end),[],2);
    TR = triangulation(Tri,ref_points);
    triplot(TR);
    pause(0.1)
end


final_points = Y(end,:);
real_points = reshape(final_points(1:2*N),[],2);

TR = triangulation(Tri,real_points);
triplot(TR);
end


function dxdt = cell_forces(t,x)
global gamma alpha s0 F N M E r_rec  t_rec T ;
dxdt = zeros(4*N,1);


t_rec = [t_rec; t];
r_rec = [r_rec; zeros(1,M)];

while t_rec(1)+T < t
    t_rec = t_rec(2:end);
    r_rec = r_rec(2:end,:);
end     %removes entries older than value T
K = size(r_rec,1);%the height of r_ij recordings
for i = 1:M;
    edge = E(i,:);
    r_ij = [x(edge(1)) x(edge(1)+N)]-[x(edge(2)) x(edge(2)+N)];
    rho_ij = [x(edge(1)+2*N) x(edge(1)+3*N)]-[x(edge(2)+2*N) x(edge(2)+3*N)];
    real_force = (norm(r_ij)-norm(rho_ij))*r_ij/norm(r_ij);
    r_rec(K,i) = norm(r_ij);%records rij
    if t_rec(end)>t_rec(1)
        r_av = trapz(t_rec,r_rec(:,i))/(t_rec(end)-t_rec(1));
    else
        r_av = r_rec(end,i);
    end
    
    ref_force = gamma*((alpha*(norm(r_ij)-s0)+(1-alpha)*(norm(rho_ij)-norm(r_av)))*rho_ij/norm(rho_ij));
    dxdt(edge(1)) = dxdt(edge(1))-real_force(1);%real xvalues
    dxdt(edge(2)) = dxdt(edge(2))+real_force(1);%real xvalues
    dxdt(edge(1)+N) = dxdt(edge(1)+N)-real_force(2);%real yvalues
    dxdt(edge(2)+N) = dxdt(edge(2)+N)+real_force(2);%real yvalues
    
    dxdt(edge(1)+2*N) = dxdt(edge(1)+2*N)-ref_force(1);%ref xvalues
    dxdt(edge(2)+2*N) = dxdt(edge(2)+2*N)+ref_force(1);%ref xvalues
    dxdt(edge(1)+3*N) = dxdt(edge(1)+3*N)-ref_force(2);%ref yvalues
    dxdt(edge(2)+3*N) = dxdt(edge(2)+3*N)+ref_force(2);%ref yvalues
end
dxdt = dxdt +F;
end
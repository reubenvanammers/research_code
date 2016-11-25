function dxdt = cell_forces(t,x)
%for spring based model, implements internal force calculation and an
%external force. 
global eta alpha s0 F N M E r_rec  t_rec T k ;
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
    real_force = k*(norm(r_ij)-norm(rho_ij))*r_ij/norm(r_ij);
    r_rec(K,i) = norm(r_ij);%records rij
    if t_rec(end)>t_rec(1)
        r_av = trapz(t_rec,r_rec(:,i))/(t_rec(end)-t_rec(1));
    else
        r_av = r_rec(end,i);
    end
    
    ref_force = k*eta*((alpha*(norm(r_ij)-s0)+(1-alpha)*(norm(rho_ij)-norm(r_av)))*rho_ij/norm(rho_ij));
    dxdt(edge(1)) = dxdt(edge(1))-real_force(1);%real xvalues
    dxdt(edge(2)) = dxdt(edge(2))+real_force(1);%real xvalues
    dxdt(edge(1)+N) = dxdt(edge(1)+N)-real_force(2);%real yvalues
    dxdt(edge(2)+N) = dxdt(edge(2)+N)+real_force(2);%real yvalues
    
    dxdt(edge(1)+2*N) = dxdt(edge(1)+2*N)-ref_force(1);%ref xvalues
    dxdt(edge(2)+2*N) = dxdt(edge(2)+2*N)+ref_force(1);%ref xvalues
    dxdt(edge(1)+3*N) = dxdt(edge(1)+3*N)-ref_force(2);%ref yvalues
    dxdt(edge(2)+3*N) = dxdt(edge(2)+3*N)+ref_force(2);%ref yvalues
end
dxdt = dxdt +F(t,x);
end
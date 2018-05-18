function dxdt = stress_force_sync(t,dxdt,external_force)
%adds external force to vertices marked in movelist. averages out restoring
%force and adds it to the external force, so that all vertices move at the
%same speed. 
global movelist restoring_rec 


[V,~] = matricize(dxdt);
N = size(V,1);
Vx = V(:,1);
restoring_forces = Vx(movelist,:);
av_restoring_force = mean(restoring_forces);
net_force = external_force+av_restoring_force;
F = columnize([net_force*movelist zeros(N,1)],zeros(N,2));
restoring_rec = [restoring_rec; abs(av_restoring_force)];

dxdt = dxdt - [movelist; movelist; zeros(2*N,1)].*dxdt;

dxdt = dxdt +F;

% av_restoring_force
% restoring_forces
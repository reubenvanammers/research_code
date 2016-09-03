function F = stress_force_sync(t,dxdt)
global movelist N external_force restoring_rec restoring_t_rec;
external_force = 0.2;
[V,~] = matricize(dxdt);
Vx = V(:,1);
restoring_forces = Vx(movelist,:);
av_restoring_force = mean(restoring_forces);
net_force = external_force+av_restoring_force;
F = columnize([net_force*movelist zeros(N,1)],zeros(N,2));
restoring_rec = [restoring_rec; abs(av_restoring_force)];
restoring_t_rec = [restoring_t_rec;t];
% av_restoring_force
% restoring_forces
function dxdt = strain_force(t,dxdt,strain,initial_length,t_strain_end)
%adds external force to vertices marked in movelist. averages out restoring
%force and adds it to the external force, so that all vertices move at the
%same speed. 
global movelist restoring_rec stress_rec
h = 0.0001;

[V,~] = matricize(dxdt);
N = size(V,1);
Vx = V(:,1);
restoring_forces = Vx(movelist,:);
av_restoring_force = mean(restoring_forces);

if t < t_strain_end
    current_derivative = initial_length*(strain(t+h)-strain(t))/h;
    F = columnize([current_derivative*movelist,zeros(N,1)],zeros(N,2));
else
    current_derivative = av_restoring_force;
    F = columnize([av_restoring_force*movelist zeros(N,1)],zeros(N,2));
end
dxdt = dxdt - [movelist; movelist; zeros(2*N,1)].*dxdt;
dxdt = dxdt +F;

restoring_rec = [restoring_rec; -av_restoring_force];
stress_rec = [stress_rec; current_derivative-av_restoring_force];

% av_restoring_force
% restoring_forces
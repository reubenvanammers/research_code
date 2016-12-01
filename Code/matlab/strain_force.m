function F = strain_force(t,dxdt,strain,initial_length)
%adds external force to vertices marked in movelist. averages out restoring
%force and adds it to the external force, so that all vertices move at the
%same speed. 
global movelist restoring_rec stress_rec
h = 0.0001;
current_derivative = (strain(t+h)-strain(t))/h;
[V,~] = matricize(dxdt);
N = size(V,1);
Vx = V(:,1);
restoring_forces = Vx(movelist,:);
av_restoring_force = mean(restoring_forces);

F = columnize([current_derivative*initial_length*movelist,zeros(N,1)],zeros(N,2));
restoring_rec = [restoring_rec; -av_restoring_force];
stress_rec = [stress_rec; current_derivative*initial_length-av_restoring_force];
% av_restoring_force
% restoring_forces
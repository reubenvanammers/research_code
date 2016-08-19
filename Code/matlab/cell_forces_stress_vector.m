function dxdt = cell_forces_stress_vector(t,x)
global gamma alpha s0 F N r_rec  t_rec T fixlist vertex_matrix_1 vertex_matrix_2 edge_matrix restoring_rec restoring_t_rec movelist;
dxdt = zeros(4*N,1);


t_rec = [t_rec; t];


while t_rec(1)+T < t
    t_rec = t_rec(2:end);
    r_rec = r_rec(2:end,:);
end     %removes entries older than value T

cell_distance = vertex_matrix_1*x-vertex_matrix_2*x;
[real_displacement, reference_displacement] = matricize(cell_distance);
real_distance = sqrt(sum(real_displacement'.^2));
reference_distance = sqrt(sum(reference_displacement'.^2));
r_rec = [r_rec; real_distance];
if t_rec(end)>t_rec(1)
    r_av = trapz(t_rec,r_rec)./(t_rec(end)-t_rec(1));
else
    r_av = r_rec(end,:);
end %calculates average distance

real_force = ((real_distance-reference_distance)./real_distance)';
ref_force = gamma.*((alpha*(reference_distance-s0)+(1-alpha)*(reference_distance-abs(r_av)))./reference_distance)';

real_force_vector = [real_force real_force].*real_displacement;
ref_force_vector = [ref_force ref_force].*reference_displacement;
total_force = columnize(real_force_vector,ref_force_vector);
dxdt = dxdt+edge_matrix*total_force;

restoring_t_rec = [restoring_t_rec; t];
restoring_rec = [restoring_rec; max(abs(dxdt([movelist; zeros(3*N,1)]==1)))];


dxdt = dxdt +F(t,x);
dxdt = dxdt - [fixlist; zeros(3*N,1)].*dxdt;
end
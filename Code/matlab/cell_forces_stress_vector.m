function dxdt = cell_forces_stress_vector(t,x)
%implements evolution of system for cell center cell centre spring based
%model to be fed into inbuilt matlab solver. Used for creep experiment.
%Done in a vectorized manner, but otherwize same as cell_forces_stress.
global eta alpha s0 N r_rec  t_rec T fixlist vertex_matrix_1 vertex_matrix_2
global edge_matrix restoring_rec movelist counter external_force k
dxdt = zeros(4*N,1);

while t_rec(end) > t;
    while t_rec(counter)+T > t;
        counter = counter-1;
    end
    t_rec = t_rec(1:end-1);
    r_rec = r_rec(1:end-1,:);
    restoring_rec = restoring_rec(1:end-1);

end
t_rec = [t_rec; t];


while t_rec(counter)+T < t
    counter = counter +1;
end     %removes entries older than value T



cell_distance = vertex_matrix_1*x-vertex_matrix_2*x;
[real_displacement, reference_displacement] = matricize(cell_distance);
real_distance = sqrt(sum(real_displacement'.^2));
reference_distance = sqrt(sum(reference_displacement'.^2));
r_rec = [r_rec; real_distance];

t_valid = t_rec(counter:end);
r_valid = r_rec(counter:end,:);

if t_valid(end)>t_valid(1)
    r_av = trapz(t_valid,r_valid)./(t_valid(end)-t_valid(1));
else
    r_av = r_valid(end,:);
end %calculates average distance

real_force = k*((real_distance-reference_distance)./real_distance)';
ref_force = k*eta.*((alpha*(reference_distance-s0)+(1-alpha)*(reference_distance-abs(r_av)))./reference_distance)';

real_force_vector = [real_force real_force].*real_displacement;
ref_force_vector = [ref_force ref_force].*reference_displacement;
total_force = columnize(real_force_vector,ref_force_vector);
dxdt = dxdt+edge_matrix*total_force;

% restoring_t_rec = [restoring_t_rec; t];
% restoring_rec = [restoring_rec; max(abs(dxdt([movelist; zeros(3*N,1)]==1)))];


dxdt = stress_force_sync(t,dxdt,external_force);

dxdt = dxdt - [fixlist; fixlist; zeros(2*N,1)].*dxdt;
end
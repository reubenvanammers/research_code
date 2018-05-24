function dxdt = cell_vertex_stress_axis_reference_stretch(t,x)
%ode describing how the system evolves for the vertex based reference model
global C N A0_vec C0_vec lambda beta gamma M alpha t_rec C_rec A_rec T 
global fixlist movelist eta restoring_rec counter included_cell external_force
global delta axis_0 fixed_vertices length_rec angle_x_proj_rec angle_y_proj_rec stretch_value

t
while t_rec(end) > t;
    while t_rec(counter)+T > t;
        counter = counter-1;
    end
    t_rec = t_rec(1:end-1);
    A_rec = A_rec(1:end-1,:);
    C_rec = C_rec(1:end-1,:);
    restoring_rec = restoring_rec(1:end-1);
    length_rec = length_rec(1:end-1,:);
    angle_x_proj_rec = angle_x_proj_rec(1:end-1,:);    
    angle_y_proj_rec = angle_y_proj_rec(1:end-1,:);
end
t_rec = [t_rec; t];


while t_rec(counter)+T < t
    counter = counter +1;
end     %removes entries older than value T


[V,V_ref] = matricize(x);
real_cell_areas = zeros(1,M);
real_cell_circumferences = zeros(1,M);
reference_cell_areas = zeros(1,M);
reference_cell_circumferences = zeros(1,M);


length_rec = [length_rec; ones(1,M)];
angle_x_proj_rec = [angle_x_proj_rec; ones(1,M)];
angle_y_proj_rec = [angle_y_proj_rec; ones(1,M)];

for i = 1:M
    real_cell_areas(i) = cell_area(i,C,V);
    real_cell_circumferences(i) = cell_circumference(i,C,V);    
    [len,direction,vertices] = cell_axes_fix(i,C,V,fixed_vertices{i}); % give stickingess?
    length_rec(end,i) = len;
    angle_x_proj_rec(end,i)  = direction(1);
    angle_y_proj_rec(end,i)  = direction(2);
    real_axis_current{i} = {len,direction,vertices};
    reference_cell_areas(i) = cell_area(i,C,V_ref);
    reference_cell_circumferences(i) = cell_circumference(i,C,V_ref);
    [len,direction,vertices] = cell_axes_fix(i,C,V_ref,fixed_vertices{i}); % give stickingess?
    reference_axis_current{i} = {stretch_value*len,direction,vertices};
end

C_rec = [C_rec;real_cell_circumferences];
A_rec = [A_rec;real_cell_areas];


t_valid = t_rec(counter:end);
A_valid = A_rec(counter:end,:);
C_valid = C_rec(counter:end,:);
l_valid = length_rec(counter:end,:);
x_valid = angle_x_proj_rec(counter:end,:);
y_valid = angle_y_proj_rec(counter:end,:);

if t_valid(end)>t_valid(1)&&T~=0
    C_av = trapz(t_valid,C_valid,1)./(t_valid(end)-t_valid(1));
    A_av = trapz(t_valid,A_valid,1)./(t_valid(end)-t_valid(1));
    l_av = trapz(t_valid,l_valid,1)./(t_valid(end)-t_valid(1));
    x_av = trapz(t_valid,x_valid,1)./(t_valid(end)-t_valid(1));
    y_av = trapz(t_valid,y_valid,1)./(t_valid(end)-t_valid(1));
    for i = 1:M
        axis_av{i} = {l_av(1),[x_av(i) y_av(i)]/norm([x_av(i) y_av(i)]),fixed_vertices{i}};
    end
else
    C_av = real_cell_circumferences;%if T=0, can't average
    A_av = real_cell_areas;
    axis_av = real_axis_current;
    
end %calculates average Area and Circumference



real_force = vertex_internal_force_calc_axis2(C,V,included_cell,lambda,beta,gamma,delta,A0_vec,real_cell_areas,C0_vec,real_cell_circumferences,reference_axis_current,real_axis_current);
follow_force = (1-alpha).*vertex_internal_force_calc_axis2(C,V_ref,included_cell,0,0,gamma,delta,A0_vec,reference_cell_areas,C0_vec,reference_cell_circumferences,axis_av,reference_axis_current);
fix_force = alpha*vertex_internal_force_calc_axis2(C,V_ref,included_cell,lambda,beta,0,0,A0_vec,reference_cell_areas,C0_vec,reference_cell_circumferences,axis_0,reference_axis_current);



dxdt = columnize(real_force,eta*(follow_force+fix_force));
dxdt = stress_force_sync(t,dxdt,external_force);

dxdt = dxdt - [fixlist; fixlist; zeros(2*N,1)].*dxdt;
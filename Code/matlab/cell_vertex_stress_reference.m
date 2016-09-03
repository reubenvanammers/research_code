function dxdt = cell_vertex_stress_reference(t,x)

global C connectivitylist F N A0_vec C0_vec lambda beta gamma M alpha t_rec C_rec A_rec T neighbouring_cells fixlist movelist eta restoring_rec restoring_t_rec

[V,V_ref] = matricize(x);
real_cell_areas = zeros(1,M);
real_cell_circumferences = zeros(1,M);
reference_cell_areas = zeros(1,M);
reference_cell_circumferences = zeros(1,M);



for i = 1:M
    real_cell_areas(i) = cell_area(i,C,V);
    real_cell_circumferences(i) = cell_circumference(i,C,V);
    reference_cell_areas(i) = cell_area(i,C,V_ref);
    reference_cell_circumferences(i) = cell_circumference(i,C,V_ref);
end

C_rec = [C_rec;real_cell_circumferences];
A_rec = [A_rec;real_cell_areas];
t_rec = [t_rec; t];

% valid_t = t_rec>(t-T)&t_rec<=t;
% valid_t = valid_t';
while t_rec(1)+T < t
    t_rec = t_rec(2:end);
    C_rec = C_rec(2:end,:);
    A_rec = A_rec(2:end,:);
end     %removes entries older than value T

% C_rec_valid = C_rec(valid_t,:);
% A_rec_valid = A_rec(valid_t,:);
% t_rec_valid = t_rec(valid_t);
if t_rec(end)>t_rec(1)
    C_av = trapz(t_rec,C_rec,1)./(t_rec(end)-t_rec(1));
    A_av = trapz(t_rec,A_rec,1)./(t_rec(end)-t_rec(1));
else
    C_av = real_cell_circumferences;%if T=0, can't average
    A_av = real_cell_areas;
end %calculates average Area and Circumference


real_force = vertex_internal_force_calc(connectivitylist,C,V,neighbouring_cells,lambda,beta,gamma,reference_cell_areas,reference_cell_circumferences);
follow_force = (1-alpha).*vertex_internal_force_calc(connectivitylist,C,V_ref,neighbouring_cells,lambda,beta,gamma,A_av,C_av);
fix_force = alpha*vertex_internal_force_calc(connectivitylist,C,V_ref,neighbouring_cells,lambda,beta,gamma,A0_vec,C0_vec);
t;



dxdt = columnize(real_force,eta*(follow_force+fix_force));
external_force = F(t,dxdt);
dxdt = dxdt - [movelist; zeros(3*N,1)].*dxdt;
dxdt = dxdt +external_force;
dxdt = dxdt - [fixlist; zeros(3*N,1)].*dxdt;
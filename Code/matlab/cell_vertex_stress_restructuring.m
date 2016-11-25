function dxdt = cell_vertex_stress_restructuring(t,x)
%ode describing how the system evolves for the vertex based reference model
%Dont use: use intermediate timesteps with other solvers
global C F N A0_vec C0_vec lambda beta gamma M alpha t_rec C_rec A_rec T 
global fixlist movelist eta restoring_rec counter included_cell cell_history external_force

[V,V_ref] = matricize(x);
%dmin = 0.2;
%dsep = dmin*1.5;

while t_rec(end) > t;
    cell_history =cell_history(1:end-1);
    C2 = cell_history{end};
    if ~ isequal(C,C2)
        included_cell = cell_inclusion(V,C2);
    end
    C=C2;
    t_rec = t_rec(1:end-1);
    A_rec = A_rec(1:end-1,:);
    C_rec = C_rec(1:end-1,:);
    restoring_rec = restoring_rec(1:end-1);
    if counter <= length(t_rec)
        while t_rec(counter)+T > t;
            counter = counter-1;
        end
    end
end
t_rec = [t_rec; t];


while t_rec(counter)+T < t
    counter = counter +1;
end     %removes entries older than value T


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



t_valid = t_rec(counter:end);
A_valid = A_rec(counter:end,:);
C_valid = C_rec(counter:end,:);


if t_valid(end)>t_valid(1)&&T~=0
    C_av = trapz(t_valid,C_valid,1)./(t_valid(end)-t_valid(1));
    A_av = trapz(t_valid,A_valid,1)./(t_valid(end)-t_valid(1));
else
    C_av = real_cell_circumferences;%if T=0, can't average
    A_av = real_cell_areas;
end %calculates average Area and Circumference


real_force = vertex_internal_force_calc(C,V,included_cell,lambda,beta,gamma,reference_cell_areas,real_cell_areas,reference_cell_circumferences,real_cell_circumferences);
follow_force = (1-alpha).*vertex_internal_force_calc(C,V_ref,included_cell,lambda,beta,gamma,A_av,reference_cell_areas,C_av,reference_cell_circumferences);
fix_force = alpha*vertex_internal_force_calc(C,V_ref,included_cell,lambda,beta,gamma,A0_vec,reference_cell_areas,C0_vec,reference_cell_circumferences);




dxdt = columnize(real_force,eta*(follow_force+fix_force));
extforce = F(t,dxdt,external_force);
dxdt = dxdt - [movelist; movelist; zeros(2*N,1)].*dxdt;

dxdt = dxdt +extforce;
dxdt = dxdt - [fixlist; fixlist; zeros(2*N,1)].*dxdt;
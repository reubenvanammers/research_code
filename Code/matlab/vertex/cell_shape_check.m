function cell_shape_check(gridsize,leftarea,rightarea)
%Test function for checking behaviour of vertex model
if nargin ==2
    leftarea = 1;
    rightarea =1;
end
global V C connectivitylist F N A0_vec C0_vec lambda beta gamma M included_cell 
sidelength = 1/sqrt(3);
A0=sqrt(27)/2*(sidelength.^2);
C0 = 6*sidelength;
lambda = 1;
beta = 1;
gamma = 0;
[V,C,connectivitylist] = hexgrid_voronoi(gridsize);
included_cell = cell_inclusion(V,C);
N= length(V);
M = length(C);

scalefactors = ones(M,1);
m = min(V(:,1));
leftscale = V(:,1) <m+0.1;
m = max(V(:,1));
rightscale =  V(:,1) >m-0.1;%plus minus 0.1 is for minor discrepancies

for i = 1:M
    for j = find(leftscale)'
        
        if connectivitylist(i,j) == 1
            scalefactors(i) = scalefactors(i)*leftarea;
        end
    end
    for j = find(rightscale)'
        if connectivitylist(i,j) == 1
            scalefactors(i) = scalefactors(i)*rightarea;
        end
    end
end
A0_vec = scalefactors.^2*A0;
C0_vec = scalefactors*C0;

V_ref = V;
V_vec = columnize(V,V_ref);
V_vec = V_vec(2*N+1:end);%ignore reference cells for now
F=0;
tend = 1000;


options = odeset('RelTol',1e-3,'AbsTol',1e-6);
[Time,Y] = ode15s(@cell_vertex_stress,0:0.2:tend,V_vec,options);
final_hex = Y(end,:)';
[V,~] = matricize([final_hex;final_hex]);

figure

for i = 1:length(C)
    patch(V(C{i},1),V(C{i},2),i,'FaceColor','w'); % draws hexagons
end
title('final hexagons')

circ = circularity(C,V);
figure
histogram(circ,10)
figure
hex_vis(Time,Y,C)
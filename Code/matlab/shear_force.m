function F = shear_force(t,x)
%adds shear force to all cells
alpha = 0.05;
T = 10;
F=0;
if t<T;
    [R,P] = matricize(x);
    realforce = R*[0 alpha; 0 0];
    refforce = P*[0 0; 0 0];
    F = columnize(realforce,refforce);
end
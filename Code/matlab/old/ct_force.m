function F = ct_force(t,x)
%constant tension compression force based on cell positions. Not really
%meaningful
alpha = 0.05;
T = 15;
F=0;
if t<T;
    [R,P] = matricize(x);
    realforce = R*[alpha 0; 0 0];
    refforce = P*[0 0; 0 0];
    F = columnize(realforce,refforce);
end
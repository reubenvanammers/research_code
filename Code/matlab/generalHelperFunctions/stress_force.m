function F=stress_force(t,x)
%adds external force to vertices marked in movelist
global movelist N external_force;
external_force = 0.2;
F = columnize([external_force*movelist zeros(N,1)],zeros(N,2));
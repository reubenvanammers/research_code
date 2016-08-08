function F=stress_force(t,x)
global movelist N;
alpha = 0.2;
F = columnize([alpha*movelist zeros(N,1)],zeros(N,2));
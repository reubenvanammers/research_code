function force = instantforce_2D(t)
global N
if t<0.1
    force = [-20;zeros(N-2,1);20];
else
    force = zeros(N,1);
end
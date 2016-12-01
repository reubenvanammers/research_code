function reference_ode()
%can probably ignore
global N s0 F alpha
N=5;
F= @instantforce
cell_positions = (0:(N-1))';
ref_cell_positions = (0:(N-1))';
alpha = 3;%scaling parameter
total_positions = [cell_positions;ref_cell_positions];%list of all cell locations, real and reference
% default_postiion = 0:(N-1);
s0 = 1; %default length
tend = 6;% how long the simulation goes for
[T,Y] = ode15s(@cell_forces,[0 tend],total_positions);

plot(T,Y(:,1:N));
xlabel('t')
title('Real Cells')
figure
plot(T,Y(:,N+1:2*N));
xlabel('t')
title('Reference Cells')

end
function dxdt = cell_forces(t,x)

global N s0 F alpha;
dxdt = (1:2*N)';
u = diff(x(1:N))-diff(x(N+1:2*N));
dxdt(1:N) =  [u;0]-[0;u] +F(t);
w = diff(x((N+1):2*N))-s0;
dxdt(N+1:2*N) = [w;0]-[0;w]-alpha*([u;0]-[0;u]);

end
%old 1d simulation without reference cells
function noreference_ode()
global N s0
N=5;
cell_positions = (2*(0:(N-1))).';
% default_postiion = 0:(N-1);
s0 = 1; %default length
tend = 50;% how long the simulation goes for
[T,Y] = ode45(@cell_forces,[0 tend],cell_positions);



plot(T,Y);
xlabel('t')


end

function dxdt = cell_forces(t,x)

global N s0

extensions = diff(x)-s0;
dxdt =  -[extensions(1);extensions] + [extensions;0];

end
    
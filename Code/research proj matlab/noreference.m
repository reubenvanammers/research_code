N=5;
cell_positions = 2*(0:(N-1));
% default_postiion = 0:(N-1);
s0 = 1; %default length
tend = 50;% how long the simulation goes for



deltat = 0.001;
for t = 1:1:(tend/(deltat));
    extensions = diff(cell_positions(t,:))-1;
    forces = - [extensions(1) extensions] + [extensions 0];
    cell_positions(t+1,:) = cell_positions(t,:) +deltat*forces;
end
plot(0:deltat:tend,cell_positions);
xlabel('t')
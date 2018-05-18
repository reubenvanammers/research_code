function [Time,Y] = cell_intermediate_restructuring(fhandle, tend, v0, varargin)
%Restructures cells when seperation greater than some amount 
global cell_history C included_cell cell_t_history
cell_history = {C};
cell_t_history = 0;
restructuring_time = 1; %How often algo checks if need to swap
h = 0.2; %delta t
Y = v0';
Time = 0;
v = v0;

dmin = 0.18;
dsep = dmin*1.01;

for i = 0:restructuring_time:(tend-restructuring_time)
    i
    [Time0,Y0] = ode45(fhandle,i:h:(i+restructuring_time),v,varargin{:});
    Time0 = Time0(2:end); Time = [Time;Time0];
    Y0 = Y0(2:end,:); Y = [Y;Y0]; 
    [V,V_ref] = matricize(Y(end,:)');
    [~,short_edges] = cell_edge_lengths(C,V,dmin);
    while size(short_edges,1) > 0
        [C,V,V_ref] = t1swap(short_edges(1,1),short_edges(1,2),dsep,C,V,V_ref);
        included_cell = cell_inclusion(V,C);
        [~,short_edges] = cell_edge_lengths(C,V,dmin);
    end
    cell_history{size(cell_history,2)+1} = C;
    cell_t_history = [cell_t_history ; i+restructuring_time];
    v = columnize(V,V_ref);
end

function [circumferences,short_edges] = cell_edge_lengths(C,V,dmin)
short_edges = [];
circumferences = zeros(1,length(C));
for l = 1:length(C)
    vertex_list = C{l};
    cellsize = length(vertex_list);
    edge_lengths = zeros(1,cellsize);
    for i = 1:cellsize;
        edge_lengths(i) = norm(V(vertex_list(mod(i,cellsize)+1),:)-V(vertex_list(i),:));
        if edge_lengths(i) < dmin
            short_edges = [short_edges; C{l}(i) C{l}(mod(i,cellsize)+1)];
        end
    end
     circumferences(l) = sum(edge_lengths);
end
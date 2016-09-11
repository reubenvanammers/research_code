function [vertex_matrix_1,vertex_matrix_2, edge_matrix] = edge_matrix_create(E,N)
%creates matrices used to convert lists of points into edges (first two
%matrices, and edge matrix is used to convert list of forces on edges to
%lists of forces on vertices. This is used in vectorized version of the
%spring model with remodelling


l = size(E,1);

v1 = zeros(l,N);
v2 = zeros(l,N);
A = zeros(l,N);
for i = 1:l
    v1(i,E(i,1)) = 1;
    v2(i,E(i,2)) = 1;
end

vertex_matrix_1 = [v1 A A A; A v1 A A; A A v1 A; A A A v1]; 
vertex_matrix_2 = [v2 A A A; A v2 A A; A A v2 A; A A A v2];

B1 = zeros(N,l);
B = zeros(N,l);
for i = 1:l
    B(E(i,1),i) = -1;
    B(E(i,2),i) = 1;
end

edge_matrix = [B B1 B1 B1; B1 B B1 B1; B1 B1 B B1; B1 B1 B1 B];
function [P,E,T] = tri_2D_Hex2()
%creates sprign mesh in a hexagonal pattern for use in stress_2d_ode, for
%creep experiments.
N=7;
c = sqrt(3)/2; %%constant
y1= (1:N')-(N+1)/2;
x1 = ((2:2:(N-1))-(1))*c;
[X1,Y1] = meshgrid(x1,y1);
A1 = X1(:);
B1 = Y1(:);

y2 = (0.5:(N+0.5))-(N+1)/2;

x2= ((1:2:N')-(1))*c;
[X2,Y2] = meshgrid(x2,y2);
A2 = X2(:);
B2 = Y2(:);

A = [A1;A2];
B = [B1;B2];

DT = delaunayTriangulation(A,B);
P = DT.Points;
T=DT.ConnectivityList;

for e  = 1:size(T,1) %replaces all triangles with [0 0 0] if they have edges of length >1.05 
    size1 = norm(P(T(e,1),:)-P(T(e,2),:));
    size2 = norm(P(T(e,3),:)-P(T(e,1),:));
    size3 = norm(P(T(e,3),:)-P(T(e,2),:));
    size4 = max([size1 size2 size3]);
    if size4 >1.05
        T(e,:) = [0 0 0];
    end
end
T( ~any(T,2), : ) = [];% removes zero entries from above loop
DT = triangulation(T,P);

E = edges(DT);
%triplot(DT);
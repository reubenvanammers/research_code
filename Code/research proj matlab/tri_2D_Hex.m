function [P,E,T] = tri_2D_Hex()
N=5;
x1= (1:N')-(N+1)/2;
y1= (1:2:N')-(N+1)/2;
[X1,Y1] = meshgrid(x1,y1);
A1 = X1(:);
B1 = Y1(:);

x2 = (0.5:(N+0.5))-(N+1)/2;
y2 = (2:2:(N-1))-(N+1)/2;
[X2,Y2] = meshgrid(x2,y2);
A2 = X2(:);
B2 = Y2(:);

A = [A1;A2];
B = [B1;B2];

DT = delaunayTriangulation(A,B);
P = DT.Points;
T=DT.ConnectivityList;
E = edges(DT);
triplot(DT);
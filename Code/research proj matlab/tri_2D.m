function [P,E,T] = tri_2D()
N=5;
x= 1:N';
y= 1:N';
[X,Y] = meshgrid(x,y);
A = X(:);
B = Y(:);
DT = delaunayTriangulation(A,B);
P = DT.Points;
T=DT.ConnectivityList;
E = edges(DT);
triplot(DT);
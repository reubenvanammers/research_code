function [P,E,T] = tri_2D_Hex2(gridsize,c2)
%creates sprign mesh in a hexagonal pattern for use in stress_2d_ode, for
%creep experiments.

if nargin <2
    c2 = 1; %%constant
end

if nargin < 1
    nx = 7;
    ny = 8;
else
    nx = gridsize(1);
    ny = gridsize(2);
end
c=  sqrt(3)/2;

if mod(nx,2)==0
    nx2 = nx+1;
else
    nx2 = nx;
end
ny2 = ny-1;

y1= (1:(ny2))'-((ny2)+1)/2;
x1 = ((2:2:((nx2)-1))-(1))*c;
[X1,Y1] = meshgrid(x1,y1);
A1 = X1(:);
B1 = Y1(:);

y2 = (0.5:((ny2)+0.5))-(ny2+1)/2;

x2= ((1:2:(nx2))'-1)*c;
[X2,Y2] = meshgrid(x2,y2);
A2 = X2(:);
B2 = Y2(:);

A = [A1;A2];
B = [B1;B2];

DT = delaunayTriangulation(A,B);

P = DT.Points*c2;
T=DT.ConnectivityList;



for e  = 1:size(T,1) %replaces all triangles with [0 0 0] if they have edges of length >1.05 
    size1 = norm(P(T(e,1),:)-P(T(e,2),:));
    size2 = norm(P(T(e,3),:)-P(T(e,1),:));
    size3 = norm(P(T(e,3),:)-P(T(e,2),:));
    size4 = max([size1 size2 size3]);
    if size4 >1.05*c2
        T(e,:) = [0 0 0];
    end
end

if mod(nx,2)==0 %removes rightmost entries for even xvalues ; 
    initial_max = max(P(:,1));
    maxlist =  find(P(:,1) ==initial_max);
    [T,P] = cell_clear(T,P,maxlist);
end

T( ~any(T,2), : ) = [];% removes zero entries from above loop
warning('off','MATLAB:triangulation:PtsNotInTriWarnId')
DT = triangulation(T,P);

E = edges(DT);
% figure
% hold on
% triplot(DT);
% for i = 1:length(A)
% %    axis off
%     plot(A(i),B(i),'r.','MarkerSize',10)
% %    circle(A(i),B(i),0.65*c)
% end


end


function [T2,P2] = cell_clear(T,P,maxlist)
k=0;
T2 = T;
for i=1:length(P)
    if any(i==maxlist)
        for n = 1:size(T2,1);
            for m = 1:3
                if T(n,m) > i% && T2(n,m) ~= 0
                    T2(n,m) = T2(n,m)-1;
                end
                if T(n,m) == i
                    T2(n,:) = [0 0 0]; 
                end
            end
        end
    else
        k = k+1;
        P2(k,:) = P(i,:); 
    end                 
end
end
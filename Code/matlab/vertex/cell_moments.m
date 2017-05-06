function [Ixx, Iyy, Ixy,Cx,Cy,A,circumference] = cell_moments(l,C,V)

vertex_list = C{l};

Ixx = 0;
Iyy = 0;
Ixy = 0;
A = 0;
Cx = 0;
Cy = 0;

circumference = 0;

nk = length(vertex_list);
for i = 1:nk;
    x = V(vertex_list(i),1);
    y = V(vertex_list(i),2);
    xp1 = V(vertex_list(mod(i,nk)+1),1);
    yp1 = V(vertex_list(mod(i,nk)+1),2);
    
    
    Cx = Cx + (x+xp1)*(x*yp1-xp1*y);
    Cy = Cy + (y+yp1)*(x*yp1-xp1*y);
    A = A + x*yp1-xp1*y;
    circumference = circumference + norm([xp1, yp1]-[x,y]);

    Ixx = Ixx + (y^2+y*yp1+yp1^2)*(x*yp1-xp1*y);
    Iyy = Iyy + (x^2+x*xp1+xp1^2)*(x*yp1-xp1*y);
    Ixy = Ixy + (x*yp1+2*x*y+2*xp1*yp1+xp1*y)*(x*yp1-xp1*y);
end
A = A/2;
Cx = Cx/(6*A);
Cy = Cy/(6*A);
Ixx = Ixx/12;
Iyy = Iyy/12;
Ixy = -Ixy/24;



Ixx-Cy*Cy*A
Iyy-Cx*Cx*A
Ixy+Cx*Cy*A
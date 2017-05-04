function [Eig1, Eig2,Evec1,Evec2] = principle_axes(Ixx, Iyy, Ixy,Cx,Cy,A)


Ixx = Ixx-Cy*Cy*A;
Iyy = Iyy-Cx*Cx*A;
Ixy = Ixy+Cx*Cy*A;

SC = Ixx+Iyy
Ixx = Ixx/SC;
Ixy = Ixy/SC;
Iyy = Iyy/SC;

E2=Ixy^2-Ixx*Iyy;
E3=1+4*E2;
Eig1 = (1+sqrt(E3))/2;
Eig2 = (1-sqrt(E3))/2;

% Eig1 = 0.5*(Iyy+Ixx+(Ixx^2+Iyy^2-4*Ixy^2-2*Ixx*Iyy)^0.5);
% Eig2 = 0.5*(Iyy+Ixx-(Ixx^2+Iyy^2-4*Ixy^2-2*Ixx*Iyy)^0.5);

Evec1 = [Iyy-Ixx-(Ixx^2+Iyy^2-4*Ixy^2-2*Ixx*Iyy)^0.5,Ixy];
Evec2 = [Iyy-Ixx+(Ixx^2+Iyy^2-4*Ixy^2-2*Ixx*Iyy)^0.5,Ixy];
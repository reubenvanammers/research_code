function [V,C,connectivitylist, neighbouring_cells] = hexgrid_voronoi()
%creates hexagonal list of vertices V, cell of cells C, and matrix of
%connetivity between the two
c = sqrt(3) / 2;
[X,Y] = meshgrid(0:1:13);
n = size(X,1);
X = c * X;
Y = Y + repmat([0 0.5],[n,n/2]);

% Plot the hexagonal mesh, including cell borders
%[XV YV] = voronoi(X(:),Y(:)); plot(XV,YV,'b-')
%axis equal, axis([10 20 10 20]), zoom on
[V,C] = voronoin([X(:),Y(:)]);
C=C(cellfun(@no_ones,C));
figure
for i = 1:length(C)
    patch(V(C{i},1),V(C{i},2),i,'FaceColor','w'); % draws hexagons
end
title('Initial Hexagons')
V(1,:) = [0 0];

for i = 1:length(C)%makes all hexagons counter clockwise
    x = V([C{i},C{i}(1)],1);
    y = V([C{i},C{i}(1)],2);
    if ispolycw(x,y)==true;
        C{i} = fliplr(C{i});
    end
end

N = length(C);%number of cells
connectivitylist = zeros(N,length(V));
for v = 2:length(V)
    for c = 1:N
        if any(C{c}==v)
            connectivitylist(c,v) = 1;
        end
    end
end

for i = 2:length(V);
    for l = 1:length(C);
        if connectivitylist(l,i) ==1;
            j1 = C{l}(circshift(C{l}==i,-1,2));
            j2 = C{l}(circshift(C{l}==i,0,2));
            j3 = C{l}(circshift(C{l}==i,1,2));
            neighbouring_cells{i}{l} = [j1 j2 j3];
        end
    end
end
        
% for v=1:length(V) %>?????? Maybe????
%     if sum(connectivitylist(:,v)) == 1
%         connectivitylist(:,v) = zeros(N,1);
%     end
% end
end
%connectivitylist(:,2:end);%removes point at infinity


function bool = no_ones(cell)
if (any(cell==1) || length(cell) <6)%removes non-hexagonal entries
    bool = false;
else
    bool = true;
end
end
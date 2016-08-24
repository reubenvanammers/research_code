Rad3Over2 = sqrt(3) / 2;
[X Y] = meshgrid(0:1:3);
n = size(X,1);
X = Rad3Over2 * X;
Y = Y + repmat([0 0.5],[n,n/2]);

% Plot the hexagonal mesh, including cell borders
%[XV YV] = voronoi(X(:),Y(:)); plot(XV,YV,'b-')
%axis equal, axis([10 20 10 20]), zoom on
[V,C] = voronoin([X(:),Y(:)]);
N = length(C)%number of cells
connectivitylist = zeros(N,length(V));
for v = 2:length(V)
    for c = 1:N
        if any(C{c}==v)
            connectivitylist(c,v) = 1;
        end
    end
end
connectivitylist = connectivitylist(:,2:end)
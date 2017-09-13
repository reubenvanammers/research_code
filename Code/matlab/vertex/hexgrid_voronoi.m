function [V,C,connectivitylist] = hexgrid_voronoi(gridsize,sidelength)
%creates hexagonal list of vertices V, cell of cells C, and matrix of
%connetivity between the two. Neighbouring cells is list of neighbouring
%cells for each
%xsize is size of cells in x direction, ysize is size of cells in y
%direction
if nargin ==0
    xsize = 7;
    ysize = 8;
else
    xsize = gridsize(1);
    ysize = gridsize(2);
end
if nargin <2
    sidelength = 1/sqrt(3);
end
xsize = xsize +2;
ysize = ysize +1;%These two lines make the x and y sizes consisten with
%what one would expect: some cells get stripped away in order for symmetry
%and this accounts for that
c = sqrt(3)/2;
[X,Y] = meshgrid(0:1:xsize,0:1:ysize);
n = size(X);
X = c * X;
if mod(n(2),2) == 1
    k = n(2)-1;
    Y(:,1:k) = Y(:,1:k) + repmat((2*c)/sqrt(3)*[0 0.5],n(1),k/2);
else
    k = n(2);
    Y(:,1:k) = Y(:,1:k) + repmat((2*c)/sqrt(3)*[0 0.5],n(1),k/2);
end
X = X(:,1:end-1);
Y = Y(:,1:end-1);


[V,C] = voronoin([X(:),Y(:)]);

C=C(cellfun(@hex_check,C));

connectivitylist = connectivity(C,V);
[C,V] = vertex_clean(C,V,connectivitylist);

ymin = min(V(:,2));
min_v = find(V(:,2)<ymin+0.01);
C=C(cellfun(@(x) ~any(ismember(x,min_v)),C));
connectivitylist = connectivity(C,V);
[C,V] = vertex_clean(C,V,connectivitylist);
connectivitylist = connectivity(C,V);

for i = 1:length(C)%makes all hexagons counter clockwise
    x = V([C{i},C{i}(1)],1);
    y = V([C{i},C{i}(1)],2);
    if ispolycw(x,y);
        C{i} = fliplr(C{i});
    end
end
V = sqrt(3)*sidelength*V;

initial_min = min(V(:,1));
V(:,1) = V(:,1)-initial_min;


%N = length(C);%number of cells



% figure
% for i = 1:length(C)
%     patch(V(C{i},1),V(C{i},2),i,'FaceColor','w'); % draws hexagons
% end
% title('Initial Hexagons')

end
%connectivitylist(:,2:end);%removes point at infinity


function bool = hex_check(cell)%removes voronoi cells less than 6 edges length, ones containing infinity
if (any(cell==1) || length(cell) <6)%removes non-hexagonal entries
    bool = false;
else
    bool = true;
end
end

function [C2,V2] = vertex_clean(C,V,connectivitylist)
k=0;
C2 = C;
for i=1:length(V)
    if sum(connectivitylist(:,i))==0
        for n = 1:length(C2);
            for m = 1:length(C2{n})
                if C{n}(m) > i
                    C2{n}(m) = C2{n}(m)-1;
                end
            end
        end
    else
        k = k+1;
        V2(k,:) = V(i,:); 
    end                 
end
end




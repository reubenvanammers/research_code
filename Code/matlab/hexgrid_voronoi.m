function [V,C,connectivitylist, neighbouring_cells] = hexgrid_voronoi()
%creates hexagonal list of vertices V, cell of cells C, and matrix of
%connetivity between the two. Neighbouring cells is list of neighbouring
%cells for each 
c = sqrt(3) / 2;
[X,Y] = meshgrid(0:1:9);
n = size(X,1);
X = c * X;
Y = Y + repmat([0 0.5],[n,n/2]);
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
    if ispolycw(x,y)==true;
        C{i} = fliplr(C{i});
    end
end



N = length(C);%number of cells
% connectivitylist = connectivity(C,V);
% [C,V] = vertex_clean(C,V);
% connectivitylist = connectivity(C,V);
% for v = 2:length(V)
%     for c = 1:N
%         if any(C{c}==v)
%             connectivitylist(c,v) = 1;
%         end
%     end
% end

for i = 1:length(V);
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

figure
for i = 1:length(C)
    patch(V(C{i},1),V(C{i},2),i,'FaceColor','w'); % draws hexagons
end
title('Initial Hexagons')

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



function connectivitylist = connectivity(C,V)
for v = 1:length(V)
    for c = 1:length(C)
        if any(C{c}==v)
            connectivitylist(c,v) = 1;
        end
    end
end
end
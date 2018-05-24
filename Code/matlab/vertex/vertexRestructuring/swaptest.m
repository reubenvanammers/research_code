[V,C] = hexgrid_voronoi();
close all
for k = 1:10
    clf
    l =randi(length(C));
    i =randi(length(C{l}));
    [C,V] = t1swap(C{l}(i),C{l}(mod(i,length(C{l}))+1),0.3,C,V);
    for j = 1:length(C)
        patch(V(C{j},1),V(C{j},2),j,'FaceColor','w'); % draws hexagons
    end
    pause(0.3)
    
end
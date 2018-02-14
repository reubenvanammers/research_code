%Draws Hexagonal cell
function hex_show(V,C)

for i = 1:length(C)
    patch(V(C{i},1),V(C{i},2),i,'FaceColor','w'); % draws hexagons
end
title('Initial Hexagons')
% 
end
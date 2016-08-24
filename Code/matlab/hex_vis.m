function hex_vis(Time,Y,C)
figure
for i = 1:length(Time);
     clf
     title(['t = ', num2str(Time(i))])
    [V,~] = matricize([Y(i,:)';Y(i,:)']);
    for j = 1:length(C)
        patch(V(C{j},1),V(C{j},2),j,'FaceColor','w'); % draws hexagons
    end
    pause(0.03)
   
end
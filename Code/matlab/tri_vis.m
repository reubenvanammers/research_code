function tri_vis(Time,Y,Tri)


for i = 1:length(Time);
    subplot(2,1,1);
    current_points = Y(i,:);
    [R,P] = matricize(current_points);%real and reference points
    TR = triangulation(Tri,R);
    triplot(TR);
    title('Real Points')
    subplot(2,1,2);
    TR = triangulation(Tri,P);
    triplot(TR);
    title('Reference Points')
    pause(0.1)
end

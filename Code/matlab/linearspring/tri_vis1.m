function tri_vis1(Time,Y,Tri)
%outputs gif file
filename = 'pictures\expfit\creeploop\realcreep';
h =figure;
for i = 1:length(Time);
%    subplot(2,1,1);
    current_points = Y(i,:);
    [R,P] = matricize(current_points);%real and reference points
    TR = triangulation(Tri,R);
    triplot(TR);
    axis([0 13 -5 5])
    axis off

%    title('Real Points')
%    title(['t = ', num2str(Time(i))])
%     subplot(2,1,2);
%     TR = triangulation(Tri,P);
%     triplot(TR);
% %    title('Reference Points')

    % Capture the plot as an image 
    frame = getframe(h); 
    im = frame2im(frame); 
    [imind,cm] = rgb2ind(im,256); 
    % Write to the GIF File 
    if i == 1 
      imwrite(imind,cm,filename,'gif', 'Loopcount',inf,'DelayTime',0.03); 
    else 
      imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',0.03); 
    end
    %imwrite(imind,cm,[filename '-' num2str(i) '.png'],'png','BackGroundColor',1)
    print('-dpng',[filename '-' num2str(i) '.png'])
    pause(0.1)
end
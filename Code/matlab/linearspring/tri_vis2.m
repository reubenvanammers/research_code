function tri_vis2(Time,Y,Tri)
%outputs gif file
filename = 'pictures\expfit\creeploop\referencecreep';
h =figure;
for i = 1:length(Time);
    clf
    subplot(2,1,1);
    current_points = Y(i,:);
    [R,P] = matricize(current_points);%real and reference points
    hold on
    TR = triangulation(Tri,R);
    plot(12.6584*ones(1,100),linspace(-5,5),'r')
    plot(0*ones(1,100),linspace(-5,5),'k--')
    plot(7.7942*ones(1,100),linspace(-5,5),'k--')
    triplot(TR);
    hold off

    axis([-3 18 -5 5])
    %annotation('line',[12.6584 -3],[12.6584 -3])
    axis off

    title('Real Cell Centres')
    pause(0.1)
%    title(['t = ', num2str(Time(i))])
    subplot(2,1,2);
    TR = triangulation(Tri,P);
    hold on
    triplot(TR);
    plot(0*ones(1,100),linspace(-5,5),'k--')
    plot(7.7942*ones(1,100),linspace(-5,5),'k--')
    hold off
    axis([-3 18 -5 5])
    axis off
    title('Reference Cell Centres')

    % Capture the plot as an image 
    frame = getframe(h); 
    im = frame2im(frame); 
    [imind,cm] = rgb2ind(im,256); 
    % Write to the GIF File 
%     if i == 1 
%       imwrite(imind,cm,filename,'gif', 'Loopcount',inf,'DelayTime',0.03); 
%     else 
%       imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',0.03); 
%     end
    %imwrite(imind,cm,[filename '-' num2str(i) '.png'],'png')
    print('-dpng',[filename '-' num2str(i) '.png'])

    pause(0.1)
end
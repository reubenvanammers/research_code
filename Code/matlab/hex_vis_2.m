function F = hex_vis_2(Time,Y,C_hist,t_rec)
%given an output from ode45/15s as [Time,Y], with real and reference data,
%and the connectivity data for the real (and thus reference) state as 
%C, displays output as a video of how system evolved
%Hex vis will show only real cells, while hex_vis_2 will show both
figure
if nargin ==3
    C = C_hist;
end
if nargin ==4
    t_rec = [t_rec;inf];
    counter = 1;
end
for i = 1:length(Time);
    if nargin ==4
        while Time(i) >= t_rec(counter);
            C = C_hist{counter};
            counter = counter +1;
        end
    end
    clf
    [V,V_ref] = matricize(Y(i,:)');
    subplot(2,1,1);
    title(['t = ', num2str(Time(i))])
    
    for j = 1:length(C)
        patch(V(C{j},1),V(C{j},2),j,'FaceColor','w'); % draws hexagons
    end
    
    subplot(2,1,2);
    %title(['t = ', num2str(Time(i))])
    
    for j = 1:length(C)
        patch(V_ref(C{j},1),V_ref(C{j},2),j,'FaceColor','w'); % draws hexagons
    end
    F(i) =getframe;
    clf;

end

movie(F,1,60);
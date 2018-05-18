function tri_vis(Time,Y,Tri)
%given an output from ode45/15s as [Time,Y], with real and reference data,
%and the triangulation data for the real (and thus reference) state as 
%Tri, displays output as a video of how system evolved

for i = 1:length(Time);
    subplot(2,1,1);
    current_points = Y(i,:);
    [R,P] = matricize(current_points);%real and reference points
    TR = triangulation(Tri,R);
    triplot(TR);
%    title('Real Points')
    title(['t = ', num2str(Time(i))])
    subplot(2,1,2);
    TR = triangulation(Tri,P);
    triplot(TR);
%    title('Reference Points')
    pause(0.1)
end

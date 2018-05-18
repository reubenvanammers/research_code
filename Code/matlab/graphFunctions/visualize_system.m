%Takes a function handle of a given simulation, and its parameters, and
%displays a movie with progression strain, cells, restoring force for
%visualisation purposes.
function [Time,Y,strain,restoring,F] = visualize_system(fhandle,varargin)

global restoring_rec t_rec
[~,C] = hexgrid_voronoi([10 10]);
[Time,Y]=fhandle(varargin{:});
l = length(t_rec)-length(restoring_rec);
t_rec = t_rec(l+1:end);
N = size(Y,2)/4;
xvalues = Y(:,1:N);
strain = (max(xvalues,[],2)-min(xvalues(1,:)))/(max(xvalues(1,:))-min(xvalues(1,:)));
restoring_temp = restoring_rec;
[t_rec,ia,~] = unique(t_rec);
restoring_temp = restoring_temp(ia);%deletes duplicate time entries for interpolation
restoring = interp1(t_rec,restoring_temp,Time);%interpolates restoring force to be same size as time vector
figure
for i = 1:length(Time);
    clf
    [V,V_ref] = matricize([Y(i,:)']);
    subplot(2,2,1);
    title(['t = ', num2str(Time(i))])
    
    for j = 1:length(C)
        patch(V(C{j},1),V(C{j},2),j,'FaceColor','w'); % draws hexagons
    end
    
    subplot(2,2,2);
    %title(['t = ', num2str(Time(i))])
    
    for j = 1:length(C)
        patch(V_ref(C{j},1),V_ref(C{j},2),j,'FaceColor','w'); % draws hexagons
    end

     
    subplot(2,2,3)
    plot(Time(1:i),strain(1:i));
    title('strain')
    axis([0 Time(end) 0.9 1.5])
    
    subplot(2,2,4)
    plot(Time(1:i),restoring(1:i));
    title('restoring')
    axis([0 Time(end) 0.8 1.1])
    F(i) = getframe(gcf);
end
figure
[h, w, p] = size(F(1).cdata);  % use 1st frame to get dimensions
set(gcf, 'position', [0 0 w+100 h+100]);
axis off
%movie(F);
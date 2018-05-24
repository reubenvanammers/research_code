%Varies strength of axial parameter, to see effect on strain plots. Only
%has a real state, with varying default axial length - this is normally
%encoded in the reference state

gamma_vec = 0:0.5:10;
target_length_vec = 1:0.1:1.6;
straincell= cell(length(gamma_vec),length(target_length_vec));
timecell =  cell(length(gamma_vec),length(target_length_vec));
for gamma_index = 1:length(gamma_vec)
    for length_index = 1:length(target_length_vec)
        [Time,strain,Y] = strain_calc(@creep_vertex_axis_stretch,100,10,gamma_vec(gamma_index),0,target_length_vec(length_index),1,0,0,1000,[10,10],0.135);
        straincell{gamma_index,length_index} = strain;
        timecell{gamma_index,length_index} = Time;
    end
end
%%

figure
hold on
colourvec = {'r','b','g','k','c','y','m','r','b','g','k','c','y','m'};
legendcell = cell(1,length(target_length_vec));
for length_index = 1:length(target_length_vec)
    strain = zeros(1,length(gamma_vec));
    for gamma_index = 1:length(gamma_vec)
        strain(gamma_index) = straincell{gamma_index,length_index}(end);
    end
    xlabel('\gamma')
    ylabel('End Strain')
    strain;
    plot(gamma_vec,strain,colourvec{length_index})
    legendcell{length_index} = ['L_0 = ',num2str(target_length_vec(length_index))];
    %title(num2str(target_length_vec(length_index)))
end
legendflex(legendcell)
%%
k = 5;
for length_index = 3:3
    j = 1;
    figure
    hold on
    for gamma_index = 1:k:length(gamma_vec) 
        plot(timecell{gamma_index,length_index},straincell{gamma_index,length_index},colourvec{ceil(gamma_index/k)});
        legendcell{j} = ['\gamma = ',num2str(gamma_vec(gamma_index))]; 
        j = j+1;
    end
    title(['Target Length = ' num2str(target_length_vec(length_index))])
    xlabel('Time')
    ylabel('Strain')
    legendflex(legendcell)
end
        
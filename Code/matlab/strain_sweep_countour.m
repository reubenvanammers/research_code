clear all;
ramptimevec = [10 15 20];
Tvec = [ 0 1 10];
etavec = logspace(-1,0,5);
alphavec = logspace(-1,0,5);
stylevec = {'-','--','-.'};
colourvec = ['r','b','k','g','y','m','c','r','y','m','c',];%need to stop reuse of colours, temp measure
error_threshold = 10^(-3.3);
figure
hold on;
[X,Y] = meshgrid(etavec,alphavec);
contours = logspace(-10,10,21);
errorcell = cell(length(Tvec));
clear h


for T_index = 1:length(Tvec);
    T = Tvec(T_index);
    load([pwd '/workspaces/strainerror' num2str(T) '_' num2str(etavec(1)) '-' num2str(etavec(end)) '_' num2str(alphavec(1)) '-' num2str(alphavec(end)) '.mat']);
    errorcell{T_index} = error_fit;
end
legendcell = {};
for T_index = 1:length(Tvec);
    for ramptime_index = 1:length(ramptimevec)
        [~,h((T_index-1)*t+ramptime_index)] = contour(X,Y,reshape(errorcell{T_index}(ramptime_index,:,:),[g,a])',[error_threshold,error_threshold],[colourvec(ramptime_index) stylevec{T_index}],'ShowText','off');
        set(gca, 'XScale', 'log', 'YScale', 'log');
        legendcell{((T_index-1)*t+ramptime_index)} =['T = ', num2str(Tvec(T_index)), ',ramptime = ' , num2str(ramptimevec(ramptime_index))];
    end
end
legend(h,legendcell)

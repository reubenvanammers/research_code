clear all;
forcevec = logspace(-1,-0.5 ,3);
Tvec = [ 0 1 10];
etavec = logspace(-1,0,10);
alphavec = logspace(-1,0,10);
stylevec = {'-','--','-.'};
colourvec = ['r','b','k','g','y','m','c','r','y','m','c',];%need to stop reuse of colours, temp measure
error_threshold = 10^(-5);
figure
hold on;
[X,Y] = meshgrid(etavec,alphavec);
contours = logspace(-10,10,21);
errorcell = cell(length(Tvec));
clear h


for T_index = 1:length(Tvec);
    T = Tvec(T_index);
    load([pwd '/workspaces/creeperror' num2str(T) '_' num2str(etavec(1)) '-' num2str(etavec(end)) '_' num2str(alphavec(1)) '-' num2str(alphavec(end)) '.mat']);
    errorcell{T_index} = error_fit;
end
legendcell = {};
for T_index = 1:length(Tvec);
    for force_index = 1:length(forcevec)
        [~,h((T_index-1)*f+force_index)] = contour(X,Y,reshape(errorcell{T_index}(force_index,:,:),[g,a])',[error_threshold,error_threshold],[colourvec(force_index) stylevec{T_index}],'ShowText','off');
        set(gca, 'XScale', 'log', 'YScale', 'log');
        xlabel('eta');
        ylabel('alpha');
        legendcell{((T_index-1)*f+force_index)} =['T = ', num2str(Tvec(T_index)), ',force = ' , num2str(forcevec(force_index))];
    end
end
legendflex(h,legendcell)
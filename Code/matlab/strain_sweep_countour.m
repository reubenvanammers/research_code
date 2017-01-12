ramptimevec = 5:5:50;
Tvec = [0 1 10];
etavec = logspace(-1,0,5);
alphavec = logspace(-1,0,5);
stylevec = ['-','--','-.'];
colourvec = ['y','m','c','r','b','g','k','y','m','c','r'];%need to stop reuse of colours, temp measure

figure
hold on;
[X,Y] = meshgrid(etavec,alphavec);
contours = logspace(-10,10,21);
h = [];

for T_index = 1:length(Tvec);
    load([pwd '\workspaces\strainerror' num2str(T) '_' num2str(etavec(1)) '-' num2str(etavec(end)) '_' num2str(alphavec(1)) '-' num2str(alphavec(end)) '.mat']);
    for ramptime_index = 1:t
        [~,h((T_index-1)*t+time_index)] = contour(X,Y,reshape(error_fit(time_index,:,:),[g,a])',[error_threshold,error_threshold],colourvec(time_index),'ShowText','off',stylevec(T_index));
        set(gca, 'XScale', 'log', 'YScale', 'log');
    end
end
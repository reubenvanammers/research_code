function strain_time(alpha,gamma,T,tend)

global N restoring_rec restoring_t_rec
[Time,Y]=stress_2d_ode(alpha,gamma,T,tend);
xvalues = Y(:,1:N);
strain = max(xvalues,[],2)/(max(xvalues(1,:))-min(xvalues(1,:)));
figure
plot(Time,strain)
xlabel('Time')
ylabel('Strain')
title(['alpha = ',num2str(alpha),' gamma = ',num2str(gamma),' T = ',num2str(T)])
set(gcf,'PaperPositionMode','auto')
saveas(gcf,[pwd '\strainfigures\' num2str(alpha) '-' num2str(gamma) '-' num2str(T) '.png']);
saveas(gcf,[pwd '\strainfigures\' num2str(alpha) '-' num2str(gamma) '-' num2str(T) '.fig']);
figure
plot(restoring_t_rec,restoring_rec)
xlabel('Time')
ylabel('Restoring Force')
title(['alpha = ',num2str(alpha),' gamma = ',num2str(gamma),' T = ',num2str(T)])
set(gcf,'PaperPositionMode','auto')
saveas(gcf,[pwd '\restoringforcefigures\' num2str(alpha) '-' num2str(gamma) '-' num2str(T) '.png']);
saveas(gcf,[pwd '\restoringforcefigures\' num2str(alpha) '-' num2str(gamma) '-' num2str(T) '.fig']);
%close all
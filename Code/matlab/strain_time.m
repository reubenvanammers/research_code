function strain_time(alpha,gamma,T)
global N external_force
[Time,Y]=stress_2d_ode(alpha,gamma,T);
xvalues = Y(:,1:N);
strain = max(xvalues,[],2)/max(xvalues(1,:));
figure
plot(Time,strain)
xlabel('Time')
ylabel('Strain')
title(['alpha = ',num2str(alpha),' gamma = ',num2str(gamma),' T = ',num2str(T)])
set(gcf,'PaperPositionMode','auto')
saveas(gcf,[pwd '\strainfigures\' num2str(alpha) '-' num2str(gamma) '-' num2str(T) '.png']);
saveas(gcf,[pwd '\strainfigures\' num2str(alpha) '-' num2str(gamma) '-' num2str(T) '.fig']);
figure
restoring_force = -1*(diff(max(xvalues,[],2))./diff(Time)-external_force)./external_force;
plot(Time(2:end),restoring_force)
xlabel('Time')
ylabel('Restoring Force')
title(['alpha = ',num2str(alpha),' gamma = ',num2str(gamma),' T = ',num2str(T)])
set(gcf,'PaperPositionMode','auto')
saveas(gcf,[pwd '\restoringforcefigures\' num2str(alpha) '-' num2str(gamma) '-' num2str(T) '.png']);
saveas(gcf,[pwd '\restoringforcefigures\' num2str(alpha) '-' num2str(gamma) '-' num2str(T) '.fig']);
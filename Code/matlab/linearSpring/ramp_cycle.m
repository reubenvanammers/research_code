%Does an example of the cyclic loading experiment using the ramp function,
%and plots stress strain curves.
strain_min = 1.2;
strain_max = 1.6;
tend = 200;
cycle_length = 20;
alphavec = logspace(-1,0,3);
etavec =logspace(-1,0,3);
a = length(alphavec);
e = length(etavec);
T = 10;

stresscell = cell(a,e);
straincell = cell(a,e);
for alpha_index = 1:a
    for eta_index = 1:e
        [alpha_index,eta_index]
        alpha = alphavec(alpha_index);
        eta = etavec(eta_index);
        strainfunc = ramp2(strain_min,strain_max,inf,cycle_length);
        [Time,Y,Tri2,stress_rec2]=strain_2d_ode(alpha,eta,T,tend,strainfunc,inf,[10 10]);
        stresscell{alpha_index,eta_index} = stress_rec2;
        N = size(Y,2)/4;
        xvalues = Y(:,1:N);
        strain = (max(xvalues,[],2)-min(xvalues(1,:)))/(max(xvalues(1,:))-min(xvalues(1,:)));
        straincell{a,e} = strain;
    end
end
save([pwd '/workspaces/rampcycle' num2str(T) '_' num2str(etavec(1)) '-' num2str(etavec(end)) '_' num2str(alphavec(1)) '-' num2str(alphavec(end)) '.mat']);

%%
load([pwd '/workspaces/rampcycle' num2str(T) '_' num2str(etavec(1)) '-' num2str(etavec(end)) '_' num2str(alphavec(1)) '-' num2str(alphavec(end)) '.mat']);

figure
for alpha_index = 1:a
    for eta_index = 1:e
        subplot(a,e,eta_index-e*(alpha_index)+a*e)
        alpha = alphavec(alpha_index);
        eta = etavec(eta_index);
        strainfunc = ramp2(strain_min,strain_max,inf,cycle_length);
        plot(Time,stresscell{alpha_index,eta_index},'k',Time,straincell{a,e},'r')
        legend(['Stress';'Strain'])
        xlabel('Time')
        title(['alpha  = ' num2str(alphavec(alpha_index)) ' eta = ' num2str(etavec(eta_index)) ])
    end
end
figure

for alpha_index = 1:a
    for eta_index = 1:e
        subplot(a,e,eta_index-e*(alpha_index)+a*e)
        alpha = alphavec(alpha_index);
        eta = etavec(eta_index);
        strainfunc = ramp2(strain_min,strain_max,inf,cycle_length);
        plot(straincell{a,e},stresscell{alpha_index,eta_index})
        xlabel('Strain')
        ylabel('Stress')
        title(['alpha  = ' num2str(alphavec(alpha_index)) ' eta = ' num2str(etavec(eta_index)) ])

    end
end
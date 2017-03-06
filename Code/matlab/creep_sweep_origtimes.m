%Varies parameters of stress_2d_ode, and plots time strain graphs in a
%subplot.
clear all

forcevec = logspace(-1.5,-1 ,3);
%Tvec = [0 20 100];
T = 0;
etavec = logspace(-1,0,9);
alphavec = logspace(-1,0,9);
tend = 200000;
f = length(forcevec);g = length(etavec);a = length(alphavec);
L = f*g*a;

straincell = cell(1,L);
timecell = cell(1,L);
flagcell = cell(1,L);
%restoringcell = cell(1,L);
parfor_progress(L);
parfor i = 1:L%index loops over alpha, then eta, then T
    counter = i-1;
    alpha_index = mod(counter,a);
    counter = (counter-alpha_index)/a;
    eta_index = mod(counter,g);
    counter = (counter-eta_index)/g;
    force_index = counter;
    alpha = alphavec(alpha_index+1);
    force = forcevec(force_index+1);
    eta = etavec(eta_index+1);%converts linear index to alpha,eta,T
    [Time, Y,~,flag] =stress_2d_ode_maxstrain(alpha,eta,T,tend,[10,10],force);
    N = size(Y,2)/4;
    xvalues = Y(:,1:N);
    strain = (max(xvalues,[],2)-min(xvalues(1,:)))/(max(xvalues(1,:))-min(xvalues(1,:)));
    straincell{i} = strain;
    timecell{i} = Time;
    flagcell{i} = flag;
    parfor_progress;
end
parfor_progress(0);

straincell2 = cell(f,g,a);
timecell2 = cell(f,g,a);

fit1 = nan*ones(f,g,a,3,3);
fit2 = nan*ones(f,g,a,3,5);
error1 = nan*ones(f,g,a,3);
error2 = nan*ones(f,g,a,3);
error_fit_L2 = nan*ones(f,g,a,3);
error_fit_inf = nan*ones(f,g,a,3);
%save([pwd '/workspaces/creeperror' num2str(T) '_' num2str(etavec(1)) '-' num2str(etavec(end)) '_' num2str(alphavec(1)) '-' num2str(alphavec(end)) '.mat']);
%%
%load([pwd '/workspaces/creeperror' num2str(T) '_' num2str(etavec(1)) '-' num2str(etavec(end)) '_' num2str(alphavec(1)) '-' num2str(alphavec(end)) '.mat']);
maxstraincell = cell(f,g,a);
timeendcell = cell(f,g,a);
straincell3 = cell(f,g,a);
timecell3 = cell(f,g,a);
for i = 1:L
    counter = i-1;
    alpha_index = mod(counter,a);
    counter = (counter-alpha_index)/a;
    eta_index = mod(counter,g);
    counter = (counter-eta_index)/g;
    force_index = counter;
    vars = {force_index+1,eta_index+1,alpha_index+1};
    straincell2{vars{:}} = straincell{i};
    timecell2{vars{:}} = timecell{i};
    maxstraincell{vars{:}} = max(straincell2{vars{:}});
    timeendcell{vars{:}} = 2*timecell2{vars{:}}(find(straincell2{vars{:}}<(1+0.9*(maxstraincell{vars{:}}-1)),1,'last'));
    if  timeendcell{vars{:}} > tend;
        ['Not enough time calculated for variables']
        vars
    end
    x_coord_scale = [0 1]; %rescales strain to unit square for comparison
    y_coord_scale = [0 1];
    strain_scaled = (straincell2{vars{:}}-1)/(maxstraincell{vars{:}}-1);
    time_scaled = timecell2{vars{:}};%/timeendcell{vars{:}};
    
    timecell3{vars{:}} = linspace(0,timeendcell{vars{:}},1001)';
    straincell3{vars{:}} = interp1(time_scaled,strain_scaled,timecell3{vars{:}});
%     timecell3{vars{:}} = linspace(0,timeendcell{vars{:}},1001)';
%     straincell3{vars{:}} = interp1(timecell2{vars{:}},straincell2{vars{:}},timecell3{vars{:}});

end
%%
unable_to_fit = [];
for guess_index = 1:3
    for i = 1:L
        counter = i-1;
        alpha_index = mod(counter,a);
        counter = (counter-alpha_index)/a;
        eta_index = mod(counter,g);
        counter = (counter-eta_index)/g;
        force_index = counter;
        vars = {force_index+1,eta_index+1,alpha_index+1,guess_index};
    %     straincell2{vars{:}} = straincell{i};
    %     timecell2{vars{:}} = timecell{i};
        if ~flagcell{i}
            try
                [fit1(vars{:},:),error1(vars{:}),fit2(vars{:},:),error2(vars{:}),error_fit_L2(vars{:}),error_fit_inf(vars{:})] = CalculateExponentialFits(timecell3{vars{1:end-1}}',straincell3{vars{1:end-1}}',1,guess_index);
            catch
                unable_to_fit = [unable_to_fit; vars];
            end
        else
            %unable_to_fit = [unable_to_fit; vars];
        end
    end
end
for guess_index = 1:3
    for force_index = 1:f
        for eta_index = 1:g
            for alpha_index = 1:a
                vars = {force_index,eta_index,alpha_index,guess_index};
                if fit2(vars{:},3) > fit2(vars{:},5)
                    temptimecoef = fit2(vars{:},3);
                    tempcoef = fit2(vars{:},2);
                    fit2(vars{:},3) = fit2(vars{:},5);
                    fit2(vars{:},5) = temptimecoef;
                    fit2(vars{:},2) = fit2(vars{:},4);
                    fit2(vars{:},4) = tempcoef;
                end
            end
        end
    end
end
%%
coef_scale_vals = nan*ones(f,g,a,3);
coef_scale_vals1 = nan*ones(f,g,a,3);
coef_scale_vals2 = nan*ones(f,g,a,3);
time_dif_1 =  nan*ones(f,g,a,3);
time_dif_2 =  nan*ones(f,g,a,3);
time_dif_3 =  nan*ones(f,g,a,3);
%Two coefficients of the exponential in the biexponential fit that are the
%same will have a  coef_scale_value of 1, while one which has only 1
%timescale, ie one of the coefficients is zero, will have coef_scale value
%of 0
for guess_index = 1:3
    for force_index = 1:f
        for eta_index = 1:g
            for alpha_index = 1:a;
                vars = {force_index,eta_index,alpha_index,guess_index};
                coef1 = fit2(vars{:},2);
                coef2 = fit2(vars{:},4);
                coef_scale_vals(vars{:}) = 4*(coef1*coef2)/((coef1+coef2)^2);
                if coef1 > coef2
                    coef_scale_vals1(vars{:}) = 1;
                    coef_scale_vals2(vars{:}) = min(coef_scale_vals(vars{:}),1);% min(1) used for floating point errors
                else
                    coef_scale_vals2(vars{:}) = 1;
                    coef_scale_vals1(vars{:}) = min(coef_scale_vals(vars{:}),1);
                end
                %Three different ways at measuring the difference between
                %the main timescales between the one exponential and the
                %two exponential fits
                time_dif_1(vars{:})  = abs(fit2(vars{:},5)-fit1(vars{:},3));
                time_dif_2(vars{:})  = abs(fit2(vars{:},5)/fit1(vars{:},3));
                time_dif_3(vars{:})  = abs(fit2(vars{:},5)-fit1(vars{:},3))/(fit2(vars{:},5)+fit1(vars{:},3));
            end
        end
    end
end
clear straincell timecell straincell2 timecell2
save([pwd '/workspaces/creeperrortimeorig' num2str(T) '_' num2str(etavec(1)) '-' num2str(etavec(end)) '_' num2str(alphavec(1)) '-' num2str(alphavec(end)) '.mat']);
%% plots strain-time graphs and exponential fits
load([pwd '/workspaces/creeperrortimeorig' num2str(T) '_' num2str(etavec(1)) '-' num2str(etavec(end)) '_' num2str(alphavec(1)) '-' num2str(alphavec(end)) '.mat']);
viewscale = 5; %Makes subplot display viewscale*viewscale for easier viewing
guess_value = 1;

viewscale = viewscale-1;
if mod(a,viewscale)==1 && mod(g,viewscale)==1 && a>1 && g>1 %reduces amount of graphs plotted so they don'f get too small: 9*9,13*13 etc creates 5*5 subplot
    g_scale = (g-1)/viewscale;
    a_scale = (a-1)/viewscale;
    alphavec_temp = alphavec(1:a_scale:a);
    etavec_temp = etavec(1:g_scale:g);
    a_temp = length(alphavec_temp);
    g_temp = length(etavec_temp);
else
    a_scale = 1;
    g_scale = 1;
    alphavec_temp = alphavec;
    etavec_temp = etavec;
    a_temp = a;
    g_temp = g;
end

for force_index = 1:f
    figure
    hold on
    for alpha_index = 1:a_temp 
        for  eta_index = 1:g_temp
            subplot(a_temp,g_temp,eta_index-g_temp*(alpha_index)+a_temp*g_temp)
            vars = {force_index,(eta_index-1)*g_scale+1,(alpha_index-1)*a_scale+1,guess_value};
            exp1 = @(x) fit1(vars{:},1) + fit1(vars{:},2)*exp(-x/fit1(vars{:},3));
            exp2 = @(x) fit2(vars{:},1) + fit2(vars{:},2)*exp(-x/fit2(vars{:},3))+ fit2(vars{:},4)*exp(-x/fit2(vars{:},5));
            exp2_1 = @(x) fit2(vars{:},2)*exp(-x/fit2(vars{:},3));
            exp2_2 = @(x) fit2(vars{:},4)*exp(-x/fit2(vars{:},5));
            %plot(timecell3{vars{:}},straincell3{vars{:}},'r',timecell3{vars{:}},exp1(timecell3{vars{:}}),'k--',timecell3{vars{:}},exp2(timecell3{vars{:}}),'b--')
            plot(timecell3{vars{1:end-1}},straincell3{vars{1:end-1}},'r',timecell3{vars{1:end-1}},exp1(timecell3{vars{1:end-1}}),'k--',timecell3{vars{1:end-1}},exp2(timecell3{vars{1:end-1}}),'b--',timecell3{vars{1:end-1}},exp2_1(timecell3{vars{1:end-1}}),'g-.',timecell3{vars{1:end-1}},exp2_2(timecell3{vars{1:end-1}}),'y-.')
            title(['alpha = ', num2str(alphavec_temp(alpha_index)), ' eta = ', num2str(etavec_temp(eta_index)), ' force = ' num2str(forcevec(force_index))]); 
            %axis([0 1 0 1]);
        end
    end
    %SaveAsPngEpsAndFig(-1,[pwd '/pictures/strainfitting/strainplot-' num2str(T) '-' num2str(timevec(time_index))]  , 7, 7/5, 9)
end

%%
% for force_index = 1:1
%     for alpha_index = 1:4:a
%         figure
%         plot(etavec,fit2(force_index,:,alpha_index,3),'k.',etavec,fit2(force_index,:,alpha_index,5),'b.','markers',14)
%         title(['time coefficients, force = ' num2str(forcevec(force_index)) ' alpha = ' num2str(alphavec(alpha_index))])
%         set(gca,'XScale','log','YScale','log')
%         xlabel('eta')
%         ylabel('Time coefficients')
%     end
% end
% 
% %%
% for force_index = 1:1
%     for eta_index = 1:4:g
%         figure
%         plot(alphavec,reshape(fit2(force_index,eta_index,:,3),[a 1]),'k.',alphavec,reshape(fit2(force_index,eta_index,:,5),[a 1]),'b.','markers',14)
%         title(['time coefficients, force = ' num2str(forcevec(force_index)) ' eta = ' num2str(etavec(eta_index))])
%         set(gca,'XScale','log','YScale','log')
%         xlabel('alpha')
%         ylabel('Time coefficients')
%     end
% end

%%
% for force_index = 1:f
%     for alpha_index = 1:4:a
%         figure
%         hold on
%         for eta_index = 1:g
%             vars = {force_index,eta_index,alpha_index,guess_value};
%             plot(etavec(eta_index),fit2(vars{:},3),'k.','markers',14,'MarkerEdgeColor',(1-coef_scale_vals1(vars{:}))*[1 1 1])
%             plot(etavec(eta_index),fit2(vars{:},5),'b.','markers',14,'MarkerEdgeColor',[1 1 1] + coef_scale_vals2(vars{:})*[-1 -1 0])
%             title(['time coefficients, force = ' num2str(forcevec(force_index)) ' alpha = ' num2str(alphavec(alpha_index))])
%             set(gca,'XScale','log','YScale','log')
%             xlabel('eta')
%             ylabel('Time coefficients')
%         end
%     end
% end
%%
for force_index = 3:3
    for alpha_index = 1:4:a
        figure
        hold on
        %yyaxis right
        plot(etavec,reshape(error1(force_index,:,alpha_index,guess_value),[g 1]),'m-')
        plot(etavec,reshape(error2(force_index,:,alpha_index,guess_value),[g 1]),'g-')
        %ylabel('L2 error')
        for eta_index = 1:g
   %         yyaxis left
            vars = {force_index,eta_index,alpha_index,guess_value};
            plot(etavec(eta_index),fit2(vars{:},3),'k.','markers',20*coef_scale_vals1(vars{:}),'MarkerEdgeColor',(1-coef_scale_vals1(vars{:}))*[1 1 1])
            plot(etavec(eta_index),fit2(vars{:},5),'b.','markers',20*coef_scale_vals2(vars{:}),'MarkerEdgeColor',[1 1 1] + coef_scale_vals2(vars{:})*[-1 -1 0])
            plot(etavec(eta_index),fit2(vars{:},2),'kx','markers',20*coef_scale_vals1(vars{:}))
            plot(etavec(eta_index),fit2(vars{:},4),'bx','markers',20*coef_scale_vals2(vars{:}))
            plot(etavec(eta_index),fit1(vars{:},2),'rx','markers',20)
            plot(etavec(eta_index),fit1(vars{:},3),'r.','markers',20)
            title(['time coefficients, force = ' num2str(forcevec(force_index)) ' alpha = ' num2str(alphavec(alpha_index))])
            set(gca,'XScale','log','YScale','log')
            xlabel('eta')
            ylabel('Time coefficients')
        end
    end
end
%%
for force_index = 3:3
    for eta_index = 1:4:g
        figure
        hold on
        %yyaxis right
        plot(alphavec,reshape(error1(force_index,eta_index,:,guess_value),[a 1]),'m-')
        plot(alphavec,reshape(error2(force_index,eta_index,:,guess_value),[a 1]),'g-')
        %ylabel('L2 error')
        for alpha_index = 1:a
   %         yyaxis left
            vars = {force_index,eta_index,alpha_index,guess_value};
            plot(alphavec(alpha_index),fit2(vars{:},3),'k.','markers',20*coef_scale_vals1(vars{:}),'MarkerEdgeColor',(1-coef_scale_vals1(vars{:}))*[1 1 1])
            plot(alphavec(alpha_index),fit2(vars{:},5),'b.','markers',20*coef_scale_vals2(vars{:}),'MarkerEdgeColor',[1 1 1] + coef_scale_vals2(vars{:})*[-1 -1 0])
            plot(alphavec(alpha_index),fit2(vars{:},2),'kx','markers',20*coef_scale_vals1(vars{:}))
            plot(alphavec(alpha_index),fit2(vars{:},4),'bx','markers',20*coef_scale_vals2(vars{:}))
            plot(alphavec(alpha_index),fit1(vars{:},2),'rx','markers',20)
            plot(alphavec(alpha_index),fit1(vars{:},3),'r.','markers',20)
            title(['time coefficients, force = ' num2str(forcevec(force_index)) ' eta = ' num2str(etavec(eta_index))])
            set(gca,'XScale','log','YScale','log')
            xlabel('alpha')
            ylabel('Time coefficients')
        end
    end
end
%%
[X,Y] = meshgrid(etavec,alphavec);
for force_index = 1:f
    figure
    surf(X,Y,reshape(fit2(force_index,:,:,guess_value,3),[g,a])')
    xlabel('eta');
    ylabel('alpha');
    title(['Timescale 1, force = ' num2str(forcevec(force_index))])
    set(gca, 'XScale', 'log', 'YScale', 'log');
    colorbar;
end

for force_index = 1:f
    figure
    surf(X,Y,reshape(fit2(force_index,:,:,guess_value,5),[g,a])')
    xlabel('eta');
    ylabel('alpha');
    title(['Timescale 2, force = ' num2str(forcevec(force_index))])
    set(gca, 'XScale', 'log', 'YScale', 'log');
    colorbar;
end

%%

% for force_index = 1:f
%     figure
%     surf(X,Y,reshape(time_dif_1(force_index,:,:,guess_value),[g,a])')
%     xlabel('eta');
%     ylabel('alpha');
%     title(['Timedif1, force = ' num2str(forcevec(force_index))])
%     set(gca, 'XScale', 'log', 'YScale', 'log');
%     colorbar;
% end

for force_index = 1:f
    figure
    surf(X,Y,reshape(time_dif_2(force_index,:,:,guess_value),[g,a])')
    xlabel('eta');
    ylabel('alpha');
    title(['Timedif2, force = ' num2str(forcevec(force_index))])
    set(gca, 'XScale', 'log', 'YScale', 'log');
    colorbar;
end

% for force_index = 1:f
%     figure
%     surf(X,Y,reshape(time_dif_3(force_index,:,:,guess_value),[g,a])')
%     xlabel('eta');
%     ylabel('alpha');
%     title(['Timedif3, force = ' num2str(forcevec(force_index))])
%     set(gca, 'XScale', 'log', 'YScale', 'log');
%     colorbar;
% end
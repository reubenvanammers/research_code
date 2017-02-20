%Varies parameters of stra/in_2d_ode, and plots time strain graphs in a
%subplot.
clear all
endstrain = 1.5;
ramptimevec = logspace(1,2,3);
%Tvec = [0 20 100];
T = 0;
etavec = logspace(-1,0,9);
alphavec = logspace(-1,0,9);
tend = 200000;


etavec_augmented = [etavec 1000]; %adds large eta (should be infinity in theory)
%in order to get lower bound of maximum stress for each alpha and ramptime
%value, to get bounds between which 90% of the value can be found to find a
%representative time. 
t = length(ramptimevec);g = length(etavec_augmented);a = length(alphavec);
L = t*g*a;
stresscell = cell(1,L);
timecell = cell(1,L);

%flagcell = cell(1,L);
%restoringcell = cell(1,L);
parfor_progress(L);
parfor i = 1:L%index loops over alpha, then eta, then T
    counter = i-1;
    alpha_index = mod(counter,a);
    counter = (counter-alpha_index)/a;
    eta_index = mod(counter,g);
    counter = (counter-eta_index)/g;
    ramptime_index = counter;
    alpha = alphavec(alpha_index+1);
    ramptime = ramptimevec(ramptime_index+1);
    eta = etavec_augmented(eta_index+1);%converts linear index to alpha,eta,T
    [~, ~,~,stress,trec,stress_index] =strain_2d_ode_ramp(alpha,eta,T,tend,ramp(endstrain,1,ramptime),ramptime,inf,[10,10]);
    %N = size(Y,2)/4;
    %xvalues = Y(:,1:N);
    %strain = (max(xvalues,[],2)-min(xvalues(1,:)))/(max(xvalues(1,:))-min(xvalues(1,:)));
    stresscell{i} = stress(stress_index+1:end);
    timecell{i} = trec(stress_index+1:end);
    [timecell{i},ia,~] = unique(timecell{i});
    stresscell{i} = stresscell{i}(ia);%deletes duplicate time entries for interpolation
    %flagcell{i} = flag;
    parfor_progress;
end
parfor_progress(0);

stresscell2 = cell(t,g,a);
timecell2 = cell(t,g,a);


save([pwd '/workspaces/strainerrortimeorig' num2str(T) '_' num2str(etavec(1)) '-' num2str(etavec(end)) '_' num2str(alphavec(1)) '-' num2str(alphavec(end)) '.mat']);
%%
load([pwd '/workspaces/strainerrortimeorig' num2str(T) '_' num2str(etavec(1)) '-' num2str(etavec(end)) '_' num2str(alphavec(1)) '-' num2str(alphavec(end)) '.mat']);
maxstresscell2 = cell(t,g,a);
endstresscell2 = cell(t,g,a);

for i = 1:L
    counter = i-1;
    alpha_index = mod(counter,a);
    counter = (counter-alpha_index)/a;
    eta_index = mod(counter,g);
    counter = (counter-eta_index)/g;
    ramptime_index = counter;
    vars = {ramptime_index+1,eta_index+1,alpha_index+1};
    stresscell2{vars{:}} = stresscell{i};
    timecell2{vars{:}} = timecell{i}-ramptimevec(ramptime_index+1);
    %rampindex = find(timecell2{vars{:}}>0,1,'first');
    maxstresscell2{vars{:}} = stresscell2{vars{:}}(1);
    stress = stresscell2{vars{:}};
    endstresscell2{vars{:}} = stress(end);
end
rampmaxstress = ones(a,t);
for ramptime_index = 1:t
    for alpha_index =1:a
        stress = cell2mat(maxstresscell2);
        rampmaxstress(alpha_index,ramptime_index) = min(stress(ramptime_index,:,alpha_index));
    end
end
%etavec = etavec(1:end-1);
g = length(etavec); L = t*g*a;
fit1 = nan*ones(t,g,a,3);
fit2 = nan*ones(t,g,a,5);
error1 = nan*ones(t,g,a);
error2 = nan*ones(t,g,a);
error_fit = nan*ones(t,g,a);
error_fit2 = nan*ones(t,g,a);
stresscell3 = cell(t,g,a);
timecell3 = cell(t,g,a);
timeendcell = cell(t,g,a);
maxstresscell = cell(t,g,a);
endstresscell = cell(t,g,a);


for i = 1:L
    counter = i-1;
    alpha_index = mod(counter,a);
    counter = (counter-alpha_index)/a;
    eta_index = mod(counter,g);
    counter = (counter-eta_index)/g;
    ramptime_index = counter;
    vars = {ramptime_index+1,eta_index+1,alpha_index+1};
    maxstresscell{vars{:}} = maxstresscell2{vars{:}};
    endstresscell{vars{:}} = endstresscell2{vars{:}};
    %timeendcell{vars{:}} = 2*timecell2{vars{:}}(find(stresscell2{vars{:}}>(endstresscell{vars{:}}+0.1*(maxstresscell{vars{:}}-endstresscell{vars{:}})),1,'last'));
    timeendcell{vars{:}} = 2*timecell2{vars{:}}(find(stresscell2{vars{:}}>(endstresscell{vars{:}}+0.1*(rampmaxstress(alpha_index+1,ramptime_index+1)-endstresscell{vars{:}})),1,'last'));

    %Chooses end time condition as twice the time it takes for stress to
    %get 90% of the way from initial to final value
    if  timeendcell{vars{:}} > tend;
        ['Not enough time calculated for variables']
        vars
    end
    
    time_scaled = timecell2{vars{:}};%/timeendcell{vars{:}};
    stress_scaled = (stresscell2{vars{:}}-endstresscell{vars{:}})/(maxstresscell{vars{:}}-endstresscell{vars{:}});
%    stress_scaled = (stresscell2{vars{:}}-endstresscell{vars{:}})/(rampmaxstress(alpha_index+1,ramptime_index+1)-endstresscell{vars{:}});
    
    
    timecell3{vars{:}} = linspace(0,timeendcell{vars{:}},1001)';
    stresscell3{vars{:}} = interp1(time_scaled,stress_scaled,timecell3{vars{:}});
end
%%
unable_to_fit = [];
for i = 1:L
    counter = i-1;
    alpha_index = mod(counter,a);
    counter = (counter-alpha_index)/a;
    eta_index = mod(counter,g);
    counter = (counter-eta_index)/g;
    time_index = counter;
    vars = {time_index+1,eta_index+1,alpha_index+1};
%     straincell2{vars{:}} = straincell{i};
%     timecell2{vars{:}} = timecell{i};
    %if ~flagcell{i}
    try
        [fit1(vars{:},:),error1(vars{:}),fit2(vars{:},:),error2(vars{:}),error_fit(vars{:}),error_fit2(vars{:})] = CalculateExponentialFits(timecell3{vars{:}}',stresscell3{vars{:}}');
    catch
        unable_to_fit = [unable_to_fit; vars];
    end
    %else
        %unable_to_fit = [unable_to_fit; vars];
    %end
end

for ramptime_index = 1:t
    for eta_index = 1:g
        for alpha_index = 1:a
            vars = {ramptime_index,eta_index,alpha_index};
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
save([pwd '/workspaces/strainerrortimeorig' num2str(T) '_' num2str(etavec(1)) '-' num2str(etavec(end)) '_' num2str(alphavec(1)) '-' num2str(alphavec(end)) '.mat']);
%% plots strain-time graphs and exponential fits
load([pwd '/workspaces/strainerrortimeorig' num2str(T) '_' num2str(etavec(1)) '-' num2str(etavec(end)) '_' num2str(alphavec(1)) '-' num2str(alphavec(end)) '.mat']);
viewscale = 5; %Makes subplot display viewscale*viewscale for easier viewing


viewscale = viewscale-1;
if mod(a,viewscale)==1 && mod(g,viewscale)==1 && a>1 && g>1 %reduces amount of graphs plotted so they don't get too small: 9*9,13*13 etc creates 5*5 subplot
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

for time_index = 1:t
    figure
    hold on
    for alpha_index = 1:a_temp 
        for  eta_index = 1:g_temp
            subplot(a_temp,g_temp,eta_index-g_temp*(alpha_index)+a_temp*g_temp)
            vars = {time_index,(eta_index-1)*g_scale+1,(alpha_index-1)*a_scale+1};
            exp1 = @(x) fit1(vars{:},1) + fit1(vars{:},2)*exp(-x/fit1(vars{:},3));
            exp2 = @(x) fit2(vars{:},1) + fit2(vars{:},2)*exp(-x/fit2(vars{:},3))+ fit2(vars{:},4)*exp(-x/fit2(vars{:},5));
            exp2_1 = @(x) fit2(vars{:},1)+fit2(vars{:},4)+fit2(vars{:},2)*exp(-x/fit2(vars{:},3));
            exp2_2 = @(x) fit2(vars{:},1) + fit2(vars{:},4)*exp(-x/fit2(vars{:},5));
            %plot(timecell3{vars{:}},stresscell3{vars{:}},'r',timecell3{vars{:}},exp1(timecell3{vars{:}}),'k--',timecell3{vars{:}},exp2(timecell3{vars{:}}),'b--')
            plot(timecell3{vars{:}},stresscell3{vars{:}},'r',timecell3{vars{:}},exp1(timecell3{vars{:}}),'k--',timecell3{vars{:}},exp2(timecell3{vars{:}}),'b--',timecell3{vars{:}},exp2_1(timecell3{vars{:}}),'g-.',timecell3{vars{:}},exp2_2(timecell3{vars{:}}),'y-.')
            title(['alpha = ', num2str(alphavec_temp(alpha_index)), ' eta = ', num2str(etavec_temp(eta_index)), ' ramptime = ' num2str(ramptimevec(time_index))]); 
            %axis([0 1 0 1]);
        end
    end
    %SaveAsPngEpsAndFig(-1,[pwd '/pictures/strainfitting/strainplot-' num2str(T) '-' num2str(timevec(time_index))]  , 7, 7/5, 9)
end

%%
% for ramptime_index = 1:1
%     for alpha_index = 1:4:a
%         figure
%         plot(etavec,fit2(ramptime_index,:,alpha_index,3),'k.',etavec,fit2(ramptime_index,:,alpha_index,5),'b.','markers',14)
%         title(['time coefficients, ramptime = ' num2str(ramptimevec(ramptime_index)) ' alpha = ' num2str(alphavec(alpha_index))])
%         set(gca,'XScale','log','YScale','log')
%         xlabel('eta')
%         ylabel('Time coefficients')
%     end
% end
% 
% %%
% for ramptime_index = 1:1
%     for eta_index = 1:4:g
%         figure
%         plot(alphavec,reshape(fit2(ramptime_index,eta_index,:,3),[a 1]),'k.',alphavec,reshape(fit2(ramptime_index,eta_index,:,5),[a 1]),'b.','markers',14)
%         title(['time coefficients, ramptime = ' num2str(ramptimevec(ramptime_index)) ' eta = ' num2str(etavec(eta_index))])
%         set(gca,'XScale','log','YScale','log')
%         xlabel('alpha')
%         ylabel('Time coefficients')
%     end
% end
%%
coef_scale_vals = nan*ones(t,g,a);
coef_scale_vals1 = nan*ones(t,g,a);
coef_scale_vals2 = nan*ones(t,g,a);
%Two coefficients of the exponential in the biexponential fit that are the
%same will have a  coef_scale_value of 1, while one which has only 1
%timescale, ie one of the coefficients is zero, will have coef_scale value
%of 0
for ramptime_index = 1:t
    for eta_index = 1:g
        for alpha_index = 1:a;
            vars = {ramptime_index,eta_index,alpha_index};
            coef1 = fit2(vars{:},2);
            coef2 = fit2(vars{:},4);
            coef_scale_vals(vars{:}) = 4*(coef1*coef2)/((coef1+coef2)^2);
            if coef1 > coef2
                coef_scale_vals1(vars{:}) = 1;
                coef_scale_vals2(vars{:}) = coef_scale_vals(vars{:});
            else
                coef_scale_vals2(vars{:}) = 1;
                coef_scale_vals1(vars{:}) = coef_scale_vals(vars{:});
            end
        end
    end
end
%%
for ramptime_index = 1:t
    for alpha_index = 1:4:a
        figure
        hold on
        for eta_index = 1:g
            vars = {ramptime_index,eta_index,alpha_index};
            plot(etavec(eta_index),fit2(ramptime_index,eta_index,alpha_index,3),'k.','markers',14,'MarkerEdgeColor',(1-coef_scale_vals1(vars{:}))*[1 1 1])
            plot(etavec(eta_index),fit2(ramptime_index,eta_index,alpha_index,5),'b.','markers',14,'MarkerEdgeColor',[1 1 1] + coef_scale_vals2(vars{:})*[-1 -1 0])
            title(['time coefficients, ramptime = ' num2str(ramptimevec(ramptime_index)) ' alpha = ' num2str(alphavec(alpha_index))])
            set(gca,'XScale','log','YScale','log')
            xlabel('eta')
            ylabel('Time coefficients')
        end
    end
end

%%
for ramptime_index = 1:t
    for eta_index = 1:4:g
        figure
        hold on
        for alpha_index = 1:a
            vars = {ramptime_index,eta_index,alpha_index};
            plot(alphavec(alpha_index),fit2(ramptime_index,eta_index,alpha_index,3),'k.','markers',20*coef_scale_vals1(vars{:}),'MarkerEdgeColor',(1-coef_scale_vals1(vars{:}))*[1 1 1])
            plot(alphavec(alpha_index),fit2(ramptime_index,eta_index,alpha_index,5),'b.','markers',20*coef_scale_vals2(vars{:}),'MarkerEdgeColor',[1 1 1] + coef_scale_vals2(vars{:})*[-1 -1 0])
            plot(alphavec(alpha_index),fit2(ramptime_index,eta_index,alpha_index,2),'kx','markers',20*coef_scale_vals1(vars{:}))
            plot(alphavec(alpha_index),fit2(ramptime_index,eta_index,alpha_index,4),'bx','markers',20*coef_scale_vals2(vars{:}))
            plot(alphavec(alpha_index),fit1(ramptime_index,eta_index,alpha_index,2),'rx','markers',20*coef_scale_vals2(vars{:}))
            plot(alphavec(alpha_index),fit1(ramptime_index,eta_index,alpha_index,3),'r.','markers',20)

            title(['time coefficients, ramptime = ' num2str(ramptimevec(ramptime_index)) ' eta = ' num2str(etavec(eta_index))])
            set(gca,'XScale','log','YScale','log')
            xlabel('alpha')
            ylabel('Time coefficients')
        end
    end
end
%%
[X,Y] = meshgrid(etavec,alphavec);

for ramptime_index = 1:t
    figure
    surf(X,Y,reshape(fit2(ramptime_index,:,:,3),[g,a])')
    xlabel('eta');
    ylabel('alpha');
    title(['Timescale 1, ramptime = ' num2str(ramptimevec(ramptime_index))])
    set(gca, 'XScale', 'log', 'YScale', 'log');
    colorbar;
end

for ramptime_index = 1:t
    figure
    surf(X,Y,reshape(fit2(ramptime_index,:,:,5),[g,a])')
    xlabel('eta');
    ylabel('alpha');
    title(['Timescale 2, ramptime = ' num2str(ramptimevec(ramptime_index))])
    set(gca, 'XScale', 'log', 'YScale', 'log');
    colorbar;
end
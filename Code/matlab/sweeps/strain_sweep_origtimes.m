%Varies parameters of stra/in_2d_ode, and plots time strain graphs in a
%subplot, with the original times being kept so thhat they can be used for
%additional analysis.
%%
vars = {ramptime_index,eta_index,alpha_index,T_value,guess_value};

k = cell(1,7);
k{1} = {2,11,11,1,1};
k{2} = {2,1,19,1,1};
k{3} = {2,1,3,1,1};
k{4} = {2,19,19,1,1};
k{5} = {2,19,5,1,1};
k{6} = {2,19,5,4,1};
k{7} = {2,19,5,6,1};


%%
clear all

rbounds = {1 2 3};
Tbounds = {0 2 5};
etabounds = {-1 0 21};
alphabounds = {-1 0 21};

rstring = arrayfun(@(x) [' t_r = 10^{' num2str(x,3) '}'], linspace(rbounds{:}),'UniformOutput',false);
astring = arrayfun(@(x) [' \alpha = 10^{' num2str(x,3) '}'], linspace(alphabounds{:}),'UniformOutput',false);
estring = arrayfun(@(x) [' \eta = 10^{' num2str(x,3) '}'], linspace(etabounds{:}),'UniformOutput',false);
Tstring = arrayfun(@(x) [' T = 10^{' num2str(x,3) '}'], linspace(Tbounds{:}),'UniformOutput',false); Tstring = {' T = 0', Tstring{:}};




endstrain = 1.5;
ramptimevec = logspace(rbounds{:});
Tvec = [0 logspace(Tbounds{:})];
etavec = logspace(etabounds{:});
alphavec = logspace(alphabounds{:});
tend = 200000;


etavec_augmented = [etavec 1000]; %adds large eta (should be infinity in theory)
%in order to get lower bound of maximum stress for each alpha and ramptime
%value, to get bounds between which 90% of the value can be found to find a
%representative time. 
t = length(ramptimevec);g_2 = length(etavec_augmented);a = length(alphavec);g = length(etavec);
L_2 = t*g_2*a*length(Tvec); L=t*g*a*length(Tvec);
stresscell = cell(1,L_2);
timecell = cell(1,L_2);

%flagcell = cell(1,L);
%restoringcell = cell(1,L);
parfor_progress(L_2);
parfor i = 1:L_2%index loops over alpha, then eta, then T
    counter = i-1;
    T_index = mod(counter,length(Tvec));
    counter = (counter-T_index)/length(Tvec);
    alpha_index = mod(counter,a);
    counter = (counter-alpha_index)/a;
    eta_index = mod(counter,g_2);
    counter = (counter-eta_index)/g_2;
    ramptime_index = counter;
    alpha = alphavec(alpha_index+1);
    ramptime = ramptimevec(ramptime_index+1);
    eta = etavec_augmented(eta_index+1);
    T = Tvec(T_index+1);
    %converts linear index to alpha,eta,T
    if eta_index == g %if calculating reference stress curve for max and min stresses for calculating endtime, uses no delay as it takes forever to calculate
        [~, ~,~,stress,trec,stress_index] =strain_2d_ode_ramp(alpha,eta,0,tend,ramp(endstrain,1,ramptime),ramptime,inf,[10,10]);
    else
        [~, ~,~,stress,trec,stress_index] =strain_2d_ode_ramp(alpha,eta,T,tend,ramp(endstrain,1,ramptime),ramptime,inf,[10,10]);
    end
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

stresscell2 = cell(t,g_2,a,length(Tvec));
timecell2 = cell(t,g_2,a,length(Tvec));


save([pwd '/workspaces/strainerrortimeorig' num2str(Tvec(1)) '_' num2str(Tvec(end)) '_' num2str(etavec(1)) '-' num2str(etavec(end)) '_' num2str(alphavec(1)) '-' num2str(alphavec(end)) '.mat']);
%%
load([pwd '/workspaces/strainerrortimeorig' num2str(Tvec(1)) '_' num2str(Tvec(end)) '_' num2str(etavec(1)) '-' num2str(etavec(end)) '_' num2str(alphavec(1)) '-' num2str(alphavec(end)) '.mat']);
t = length(ramptimevec);g_2 = length(etavec_augmented);a = length(alphavec);g = length(etavec);
L_2 = t*g_2*a*length(Tvec); L=t*g*a*length(Tvec);
maxstresscell2 = cell(t,g_2,a,length(Tvec));
endstresscell2 = cell(t,g_2,a,length(Tvec));

for i = 1:L_2
    counter = i-1;
    T_index = mod(counter,length(Tvec));
    counter = (counter-T_index)/length(Tvec);
    alpha_index = mod(counter,a);
    counter = (counter-alpha_index)/a;
    eta_index = mod(counter,g_2);
    counter = (counter-eta_index)/g_2;
    ramptime_index = counter;
    vars = {ramptime_index+1,eta_index+1,alpha_index+1,T_index+1};
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
        rampmaxstress(alpha_index,ramptime_index) = min(stress(ramptime_index,:,alpha_index,1));
    end
end
%etavec = etavec(1:end-1);
g = length(etavec); L = t*g*a*length(Tvec);
fit1 = nan*ones(t,g,a,length(Tvec),3,3);
fit2 = nan*ones(t,g,a,length(Tvec),3,5);
error1 = nan*ones(t,g,a,length(Tvec),3);
error2 = nan*ones(t,g,a,length(Tvec),3);
error_fit_L2 = nan*ones(t,g,a,length(Tvec),3);
error_fit_inf = nan*ones(t,g,a,length(Tvec),3);
stresscell3 = cell(t,g,a,length(Tvec));
timecell3 = cell(t,g,a,length(Tvec));
timeendcell = cell(t,g,a,length(Tvec));
maxstresscell = cell(t,g,a,length(Tvec));
endstresscell = cell(t,g,a,length(Tvec));


for i = 1:L
    counter = i-1;
    T_index = mod(counter,length(Tvec));
    counter = (counter-T_index)/length(Tvec);
    alpha_index = mod(counter,a);
    counter = (counter-alpha_index)/a;
    eta_index = mod(counter,g);
    counter = (counter-eta_index)/g;
    ramptime_index = counter;
    vars = {ramptime_index+1,eta_index+1,alpha_index+1,T_index+1};
    maxstresscell{vars{:}} = maxstresscell2{vars{:}};
    endstresscell{vars{:}} = endstresscell2{vars{:}};
    timeendcell{vars{:}} = 2*timecell2{vars{:}}(find(stresscell2{vars{:}}>(endstresscell{vars{:}}+0.1*(maxstresscell{vars{:}}-endstresscell{vars{:}})),1,'last'));
    %timeendcell{vars{:}} = 2*timecell2{vars{:}}(find(stresscell2{vars{:}}>(endstresscell{vars{:}}+0.1*(rampmaxstress(alpha_index+1,ramptime_index+1)-endstresscell{vars{:}})),1,'last'));

    %Chooses end time condition as twice the time it takes for stress to
    %get 90% of the way from initial to final value
    if  timeendcell{vars{:}} > tend;
        ['Not enough time calculated for variables']
        vars
    end
    
    time_scaled = timecell2{vars{:}};%/timeendcell{vars{:}};
%    stress_scaled = (stresscell2{vars{:}}-endstresscell{vars{:}})/(maxstresscell{vars{:}}-endstresscell{vars{:}});
    stress_scaled = (stresscell2{vars{:}}-endstresscell{vars{:}})/(rampmaxstress(alpha_index+1,ramptime_index+1)-endstresscell{vars{:}});
    
    
    timecell3{vars{:}} = linspace(0,timeendcell{vars{:}},1001)';
    stresscell3{vars{:}} = interp1(time_scaled,stress_scaled,timecell3{vars{:}});
end
%%
unable_to_fit = [];
for guess_index = 1:3 
    for i = 1:L
        counter = i-1;
        T_index = mod(counter,length(Tvec));
        counter = (counter-T_index)/length(Tvec);
        alpha_index = mod(counter,a);
        counter = (counter-alpha_index)/a;
        eta_index = mod(counter,g);
        counter = (counter-eta_index)/g;
        ramptime_index = counter;
        vars = {ramptime_index+1,eta_index+1,alpha_index+1,T_index+1,guess_index};
        
    %     straincell2{vars{:}} = straincell{i};
    %     timecell2{vars{:}} = timecell{i};
        %if ~flagcell{i}
        try
            [fit1(vars{:},:),error1(vars{:}),fit2(vars{:},:),error2(vars{:}),error_fit_L2(vars{:}),error_fit_inf(vars{:})] = CalculateExponentialFits(timecell3{vars{1:end-1}}',stresscell3{vars{1:end-1}}',-1,guess_index);
        catch
            unable_to_fit = [unable_to_fit; vars];
        end
        %else
            %unable_to_fit = [unable_to_fit; vars];
        %end
    end
end
for guess_index = 1:3
    for ramptime_index = 1:t
        for eta_index = 1:g
            for alpha_index = 1:a
                for T_index = 1:length(Tvec)
                    vars = {ramptime_index,eta_index,alpha_index,T_index,guess_index};
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
end
%%
coef_scale_vals = nan*ones(t,g,a,length(Tvec),3);
coef_scale_vals1 = nan*ones(t,g,a,length(Tvec),3);
coef_scale_vals2 = nan*ones(t,g,a,length(Tvec),3);
time_dif_1 =  nan*ones(t,g,a,length(Tvec),3);
time_dif_2 =  nan*ones(t,g,a,length(Tvec),3);
time_dif_3 =  nan*ones(t,g,a,length(Tvec),3);
%Two coefficients of the exponential in the biexponential fit that are the
%same will have a  coef_scale_value of 1, while one which has only 1
%timescale, ie one of the coefficients is zero, will have coef_scale value
%of 0
for guess_index = 1:3
    for ramptime_index = 1:t
        for eta_index = 1:g
            for alpha_index = 1:a;
                for T_index = 1:length(Tvec)
                    vars = {ramptime_index,eta_index,alpha_index,T_index,guess_index};
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
end
guess_value = 1;
T_value = 1;
viewscale = 11; %Makes subplot display viewscale*viewscale for easier viewing

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
save([pwd '/workspaces/strainerrortimeorig' num2str(Tvec(1)) '_' num2str(Tvec(end)) '_' num2str(etavec(1)) '-' num2str(etavec(end)) '_' num2str(alphavec(1)) '-' num2str(alphavec(end)) '.mat']);
%% plots strain-time graphs and exponential fits
load([pwd '/workspaces/strainerrortimeorig' num2str(Tvec(1)) '_' num2str(Tvec(end)) '_' num2str(etavec(1)) '-' num2str(etavec(end)) '_' num2str(alphavec(1)) '-' num2str(alphavec(end)) '.mat']);


%%
for ramptime_index = 1:t
    figure
    hold on
    for alpha_index = 1:a_temp 
        for  eta_index = 1:g_temp
            subplot(a_temp,g_temp,eta_index-g_temp*(alpha_index)+a_temp*g_temp)
            vars = {ramptime_index,(eta_index-1)*g_scale+1,(alpha_index-1)*a_scale+1,T_value,guess_value};
            exp1 = @(x) fit1(vars{:},1) + fit1(vars{:},2)*exp(-x/fit1(vars{:},3));
            exp2 = @(x) fit2(vars{:},1) + fit2(vars{:},2)*exp(-x/fit2(vars{:},3))+ fit2(vars{:},4)*exp(-x/fit2(vars{:},5));
            exp2_1 = @(x) fit2(vars{:},1)+fit2(vars{:},4)+fit2(vars{:},2)*exp(-x/fit2(vars{:},3));
            exp2_2 = @(x) fit2(vars{:},1) + fit2(vars{:},4)*exp(-x/fit2(vars{:},5));
            %plot(timecell3{vars{:}},stresscell3{vars{:}},'r',timecell3{vars{:}},exp1(timecell3{vars{:}}),'k--',timecell3{vars{:}},exp2(timecell3{vars{:}}),'b--')
            plot(timecell3{vars{1:end-1}},stresscell3{vars{1:end-1}},'r',timecell3{vars{1:end-1}},exp1(timecell3{vars{1:end-1}}),'k--',timecell3{vars{1:end-1}},exp2(timecell3{vars{1:end-1}}),'b--',timecell3{vars{1:end-1}},exp2_1(timecell3{vars{1:end-1}}),'g-.',timecell3{vars{1:end-1}},exp2_2(timecell3{vars{1:end-1}}),'y-.')
            title([astring{(alpha_index-1)*a_scale+1}, estring{(eta_index-1)*g_scale+1}, rstring{ramptime_index}, Tstring{T_value} ]); 
            legend('Data','Single Exp', 'Two Exp', 'Short Time', 'Long Time') 
            %axis([0 1 0 1]);
        end
    end
    %SaveAsPngEpsAndFig(-1,[pwd '/pictures/strainfitting/strainplot-' num2str(T) '-' num2str(timevec(ramptime_index))]  , 7, 7/5, 9)
end

%%
%for ramptime_index = 1:t
%    for alpha_index = 1:a_temp 
%        for  eta_index = 1:g_temp
        for i = 1:length(k)
            vars = k{i}
            ramptime_index = vars{1};
            eta_index = vars{2};
            alpha_index = vars{3};
            T_index = vars{4};
            guess_index = vars{5};
            figure
            hold on
            %subplot(a_temp,g_temp,eta_index-g_temp*(alpha_index)+a_temp*g_temp)
            %vars = {ramptime_index,(eta_index-1)*g_scale+1,(alpha_index-1)*a_scale+1,T_value,guess_value};
            exp1 = @(x) fit1(vars{:},1) + fit1(vars{:},2)*exp(-x/fit1(vars{:},3));
            exp2 = @(x) fit2(vars{:},1) + fit2(vars{:},2)*exp(-x/fit2(vars{:},3))+ fit2(vars{:},4)*exp(-x/fit2(vars{:},5));
            exp2_1 = @(x) fit2(vars{:},1)+fit2(vars{:},4)+fit2(vars{:},2)*exp(-x/fit2(vars{:},3));
            exp2_2 = @(x) fit2(vars{:},1) + fit2(vars{:},4)*exp(-x/fit2(vars{:},5));
            time = timecell3{vars{1:end-1}}; stress = stresscell3{vars{1:end-1}};
            plot(time(1:10:end),stress(1:10:end),'k.',time,exp1(time),'r-',time,exp2(time),'b--','MarkerSize',4)
%             plot(timecell3{vars{1:end-1}},stresscell3{vars{1:end-1}},'r',timecell3{vars{1:end-1}},exp1(timecell3{vars{1:end-1}}),'k--',timecell3{vars{1:end-1}},exp2(timecell3{vars{1:end-1}}),'b--',timecell3{vars{1:end-1}},exp2_1(timecell3{vars{1:end-1}}),'g-.',timecell3{vars{1:end-1}},exp2_2(timecell3{vars{1:end-1}}),'y-.')
            
            title([astring{alpha_index}, estring{eta_index}, rstring{ramptime_index}, Tstring{T_index} ]); 

            %title([astring{(alpha_index-1)*a_scale+1}, estring{(eta_index-1)*g_scale+1}, rstring{ramptime_index}, Tstring{T_value} ]); 
            %legend('Data','Single Exp', 'Two Exp', 'Short Time', 'Long Time') 
            %axis([0 1 0 1]);
            xlabel('Time')
            ylabel('Stress')
            SaveAsPngEpsAndFig([pwd '/pictures/expfit/relaxation/timestress/' num2str(Tvec(T_index)) '-' num2str(ramptimevec(ramptime_index)) '-' num2str(etavec(eta_index)) '-' num2str(alphavec(alpha_index))])
            %pause
            close all
        end
%    end
%end
%%
%%
for ramptime_index = 1:t
    figure
    [X,Y] = meshgrid(etavec,alphavec);
    hold on;
    surf(X,Y,reshape(cell2mat(timeendcell(ramptime_index,:,:,T_value,guess_value)),[g,a])');
    shading interp; 
%     alpha(0.5);
     colorbar;
    %contour(X,Y,reshape(cell2mat(timeendcell(force_index,:,:)),[g,a])',contours,'ShowText','on');

    set(gca, 'XScale', 'log', 'YScale', 'log','ZScale','log');
    xlabel('\eta');
    ylabel('\alpha');
    title(['equilibriation times,' rstring{ramptime_index}]);
    SaveAsPngEpsAndFig(-1,[pwd '/pictures/expfit/relaxation/timeendsurface/' num2str(Tvec(T_value)) '-' num2str(ramptimevec(ramptime_index))]  , 7, 7/5, 9)
    %SaveAsPngEpsAndFig(-1,[pwd 'asdf']  , 7, 7/5, 9)

end
%%
% for ramptime_index = 1:t
%     for alpha_index = 1:4:a
%         figure
%         hold on
%         for eta_index = 1:g
%             vars = {ramptime_index,eta_index,alpha_index,T_index,guess_value};
%             plot(etavec(eta_index),fit2(vars{:},3),'k.','markers',14,'MarkerEdgeColor',(1-coef_scale_vals1(vars{:}))*[1 1 1])
%             plot(etavec(eta_index),fit2(vars{:},5),'b.','markers',14,'MarkerEdgeColor',[1 1 1] + coef_scale_vals2(vars{:})*[-1 -1 0])
%             title(['time coefficients, ramptime = ' num2str(ramptimevec(ramptime_index)) ' alpha = ' num2str(alphavec(alpha_index))])
%             set(gca,'XScale','log','YScale','log')
%             xlabel('eta')
%             ylabel('Time coefficients')
%         end
%     end
% end

%%
for ramptime_index = 1:t %creates line plots in order to view differing coefficient values for eta and alpha
    for alpha_index = 1:a
        figure
        hold on
%         yyaxis right
%         h(1) = plot(etavec,reshape(error1(ramptime_index,:,(alpha_index-1)*a_scale+1,T_value,guess_value),[g 1]),'m-');
%         h(2) = plot(etavec,reshape(error2(ramptime_index,:,(alpha_index-1)*a_scale+1,T_value,guess_value),[g 1]),'g-');
%         ylabel('L2 error')
%         set(gca,'XScale','log','YScale','log')

        for eta_index = 1:g
%             yyaxis left
            vars = {ramptime_index,eta_index,alpha_index,T_value,guess_value};
            h(1) = plot(etavec(eta_index),fit2(vars{:},3),'k.','markers',10*coef_scale_vals1(vars{:}));%,'MarkerEdgeColor',(1-coef_scale_vals1(vars{:}))*[1 1 1])
            h(2) = plot(etavec(eta_index),fit2(vars{:},5),'b.','markers',10*coef_scale_vals2(vars{:}));%z`,'MarkerEdgeColor',[1 1 1] + coef_scale_vals2(vars{:})*[-1 -1 0])
            h(3) = plot(etavec(eta_index),fit1(vars{:},3),'r.','markers',10);%*coef_scale_vals2(vars{:}));
%             h(6) = plot(etavec(eta_index),fit2(vars{:},2),'kx','markers',5);%*coef_scale_vals1(vars{:}))
%             h(7) = plot(etavec(eta_index),fit2(vars{:},4),'bx','markers',5);%*coef_scale_vals2(vars{:}));
%             h(8) = plot(etavec(eta_index),fit1(vars{:},2),'rx','markers',5);

            title(['time coefficients' rstring{ramptime_index} astring{(alpha_index)}])
            set(gca,'XScale','log','YScale','log')

        end
            xlabel('\eta')
            ylabel('Time coefficients')
            %SaveAsPngEpsAndFig(-1,[pwd '/pictures/expfit/relaxation/lineplot/' num2str(Tvec(T_value)) '-' num2str(ramptimevec(ramptime_index)) '-alpha=' num2str(alphavec((alpha_index-1)*a_scale+1))]  , 7, 7/5, 9)
            legend(h)
    end
end
%%
for ramptime_index = 1:t %creates line plots in order to view differing coefficient values for eta and alpha
    for eta_index = 1:g_temp
        figure
        hold on
%         yyaxis right
%         h(1) = plot(alphavec,reshape(error1(ramptime_index,(eta_index-1)*g_scale+1,:,T_value,guess_value),[a 1]),'m-');
%         h(2) = plot(alphavec,reshape(error2(ramptime_index,(eta_index-1)*g_scale+1,:,T_value,guess_value),[a 1]),'g-');
%         ylabel('L2 error')
%         set(gca,'XScale','log','YScale','log')

        for alpha_index = 1:a
%             yyaxis left
            vars = {ramptime_index,(eta_index-1)*g_scale+1,alpha_index,T_value,guess_value};
            h(3) = plot(alphavec(alpha_index),fit2(vars{:},3),'k.','markers',10*coef_scale_vals1(vars{:}));%,'MarkerEdgeColor',(1-coef_scale_vals1(vars{:}))*[1 1 1])
            h(4) = plot(alphavec(alpha_index),fit2(vars{:},5),'b.','markers',10*coef_scale_vals2(vars{:}));%z`,'MarkerEdgeColor',[1 1 1] + coef_scale_vals2(vars{:})*[-1 -1 0])
            h(5) = plot(alphavec(alpha_index),fit1(vars{:},3),'r.','markers',10);%*coef_scale_vals2(vars{:}));
%             h(6) = plot(alphavec(alpha_index),fit2(vars{:},2),'kx','markers',5);%*coef_scale_vals1(vars{:}))
%             h(7) = plot(alphavec(alpha_index),fit2(vars{:},4),'bx','markers',5);%*coef_scale_vals2(vars{:}));
%             h(8) = plot(alphavec(alpha_index),fit1(vars{:},2),'rx','markers',5);

            title(['time coefficients,' rstring{ramptime_index} estring{(eta_index-1)*g_scale+1}])
            set(gca,'XScale','log','YScale','log')

        end
            xlabel('\alpha')
            ylabel('Time coefficients')
            SaveAsPngEpsAndFig(-1,[pwd '/pictures/expfit/relaxation/lineplot/' num2str(Tvec(T_value)) '-' num2str(ramptimevec(ramptime_index)) '-eta=' num2str(etavec((eta_index-1)*g_scale+1))]  , 7, 7/5, 9)
            %legend(h)
    end
end
%%
[X,Y] = meshgrid(etavec,alphavec);
for ramptime_index = 1:t
    figure
    surf(X,Y,reshape(fit2(ramptime_index,:,:,T_value,guess_value,3),[g,a])')
    xlabel('\eta');
    ylabel('\alpha');
    title(['Timescale 1,' rstring{ramptime_index}])
    set(gca, 'XScale', 'log', 'YScale', 'log');
    colorbar;
end

for ramptime_index = 1:t
    figure
    surf(X,Y,reshape(fit2(ramptime_index,:,:,T_value,guess_value,5),[g,a])')
    xlabel('\eta');
    ylabel('\alpha');
    title(['Timescale 2,' rstring{ramptime_index}])
    set(gca, 'XScale', 'log', 'YScale', 'log');
    colorbar;
end

%%
[X,Y] = meshgrid(etavec,alphavec);

% for ramptime_index = 1:t
%     figure
%     surf(X,Y,reshape(time_dif_1(ramptime_index,:,:,T_value,guess_value),[g,a])')
%     xlabel('eta');
%     ylabel('alpha');
%     title(['Timedif1, ramptime = ' num2str(ramptimevec(ramptime_index))])
%     set(gca, 'XScale', 'log', 'YScale', 'log');
%     colorbar;
% end

for ramptime_index = 1:t
    figure
    surf(X,Y,reshape(time_dif_2(ramptime_index,:,:,T_value,guess_value),[g,a])')
    xlabel('\eta');
    ylabel('\alpha');
    title(['Timedif2,' rstring{ramptime_index}])
    set(gca, 'XScale', 'log', 'YScale', 'log');
    colorbar;
end

% for ramptime_index = 1:t
%     figure
%     surf(X,Y,reshape(time_dif_3(ramptime_index,:,:,T_value,guess_value),[g,a])')
%     xlabel('eta');
%     ylabel('alpha');
%     title(['Timedif3, ramptime = /' num2str(ramptimevec(ramptime_index))])
%     set(gca, 'XScale', 'log', 'YScale', 'log');
%     colorbar;
% end

%%
[X,Y] = meshgrid(etavec,alphavec);

for ramptime_index = 1:t
    figure
    surf(X,Y,reshape(error1(ramptime_index,:,:,T_value,guess_value),[g,a])')
    
    xlabel('\eta');
    ylabel('\alpha');
    title(['one exponential L2 error,' rstring{ramptime_index}])
    set(gca, 'XScale', 'log', 'YScale', 'log');
end


%%
[X,Y] = meshgrid(etavec,alphavec);
time_dif_2(:,:,a,:,:) = nan*time_dif_2(:,:,a,:,:);

contours = logspace(0,5,61);
for ramptime_index = 1:t %view the difference between coefficients for one and two exp fits as a surface plot
    for T_index = 1:1%length(Tvec)
        figure
        hold on
        surf(X,Y,reshape(time_dif_2(ramptime_index,:,:,T_index,guess_value),[g,a])')
        clearvars alpha
        %alpha(0.5)
        %contour(X,Y,reshape(time_dif_2(ramptime_index,:,:,T_index,guess_value),[g,a])',contours,'ShowText','on')

        xlabel('\eta');
        ylabel('\alpha');
        title(['CoefRatio,'  rstring{ramptime_index} Tstring{T_index} ])
        set(gca, 'XScale', 'log', 'YScale', 'log');

        colorbar;
        caxis([1 5])
        %SaveAsPngEpsAndFig(-1,[pwd '/pictures/expfit/relaxation/timedifsurface/' num2str(Tvec(T_index)) '-' num2str(ramptimevec(ramptime_index))]  , 7, 7/5, 9)
        %close all
    end
end


%%
stylevec = {'-','--','-.',':'};
colourvec = {[1 0 0],[0 0 1],[0 0 0],[0 1 0],[1 1 0],[1 0 1],[0 1 1],[1 0 0],[0 0 1],[0 0 0],[0 1 0]};%need to stop reuse of colours, temp measure
error_threshold = 2;
figure
hold on;
[X,Y] = meshgrid(etavec,alphavec);

two_exp_status = time_dif_2 > error_threshold;

for T_index = 1:length(Tvec)%view compound contour plot - generally too messy to use
    for ramptime_index = 1:t
        [~,h((T_index-1)*t+ramptime_index)] = contour(X,Y,reshape(time_dif_2(ramptime_index,:,:,T_index,guess_value),[g,a])',[error_threshold,error_threshold],[stylevec{T_index}],'ShowText','off','LineColor',colourvec{ramptime_index});
        set(gca, 'XScale', 'log', 'YScale', 'log');
        legendcell{((T_index-1)*t+ramptime_index)} =['T = ', num2str(Tvec(T_index)), ',ramptime = ' , num2str(ramptimevec(ramptime_index))];
        xlabel('\eta');
        ylabel('\alpha');
        title(['overall relaxation contour, Threshold = ' num2str(error_threshold)])

    end
end

legendflex(h,legendcell)

%%
clear legendcell h
colourvec = {'r','b','g'};
sizevec = {4,5,6};
for T_index = 1
    figure
    hold on
    clear h legendcell
        %view regions of two exponential bejaviour for varying ramptimes
    for ramptime_index = t:-1:1
        [~,h(ramptime_index)] = contour(X,Y,reshape(time_dif_2(ramptime_index,:,:,T_index,guess_value),[g,a])',[error_threshold,error_threshold],[stylevec{1}],'ShowText','off','LineColor',colourvec{ramptime_index});
        set(gca, 'XScale', 'log', 'YScale', 'log');
        legendcell{(ramptime_index)} =['t_r = ' , num2str(ramptimevec(ramptime_index))];
        xlabel('\eta');
        ylabel('\alpha');
        %title([%'overall relaxation contour, 
        %    'Threshold = ' num2str(error_threshold), Tstring{T_index}])
        Z = reshape(two_exp_status(ramptime_index,:,:,T_index,guess_value),[g,a])';
        plot(Z.*X,Z.*Y,['.'],'markers',sizevec{ramptime_index},'Color',colourvec{ramptime_index})
    end
    legendflex(fliplr(h),fliplr(rstring))
    
    %SaveAsPngEpsAndFig([pwd '/pictures/expfit/relaxation/biexpcontour/' num2str(Tvec(T_index))])
    %pause
    %SaveAsPngEpsAndFig(-1,[pwd '/pictures/expfit/relaxation/biexpcontour/' num2str(Tvec(T_index))]  , 7, 7/5, 9)
end


%%
clear legendcell h

colourvec = {'r','k','k','b','k','g'};
sizevec = {4,1,1,5,1,6};
for ramptime_index = 2:2
    figure
    hold on
    clear h legendcell
        %view regions of two exponential bejaviour for varying T_mem
        %averaging values
    for T_index = [6 4 1]
        [~,h(T_index)] = contour(X,Y,reshape(time_dif_2(ramptime_index,:,:,T_index,guess_value),[g,a])',[error_threshold,error_threshold],[stylevec{1}],'ShowText','off','LineColor',colourvec{T_index});
        set(gca, 'XScale', 'log', 'YScale', 'log');
        legendcell{(T_index)} =['T = ' , num2str(Tvec(T_index))];
        xlabel('\eta');
        ylabel('\alpha');
        %title(['Threshold = ' num2str(error_threshold), rstring{ramptime_index}])
        Z = reshape(two_exp_status(ramptime_index,:,:,T_index,guess_value),[g,a])';
        plot(Z.*X,Z.*Y,['.'],'markers',sizevec{T_index},'Color',colourvec{T_index} )
        %surf(Z.*X,Z.*Y,+Z)
        %colormap(colourvec{T_index})
        %alpha(0.3)
        %shading interp
    end
    
    h2 = [h(1) h(4) h(6)];
    Tstring2 ={Tstring{1} Tstring{4} Tstring{6}};
    legendflex(h2,Tstring2)
    
    %SaveAsPngEpsAndFig([pwd '/pictures/expfit/relaxation/biexpcontour2/' num2str(ramptimevec(ramptime_index))])

end


%%
%% Plots error values of one and two exponential fits
for T_index = 1
   for ramptime_index = 1:t
       for alpha_index = 1:a
           for  eta_index = 1:g
                figure
                hold on

                time = timecell3{vars{:}}; stress = stresscell3{vars{:}};
                plot(timecell3{vars{:}},stress-exp1(timecell3{vars{:}}),'r-',timecell3{vars{:}},stress-exp2(timecell3{vars{:}}),'b-',timecell3{vars{:}},zeros(size(timecell3{vars{:}})),'k--')%,timecell3{vars{:}},stress-exp3(timecell3{vars{:}}),'g-.','MarkerSize',4)            %plot(timecell3{vars{1:end-1}},stresscell3{vars{1:end-1}},'r',timecell3{vars{1:end-1}},exp1(timecell3{vars{1:end-1}}),'k--',timecell3{vars{1:end-1}},exp2(timecell3{vars{1:end-1}}),'b--',timecell3{vars{1:end-1}},exp2_1(timecell3{vars{1:end-1}}),'g-.',timecell3{vars{1:end-1}},exp2_2(timecell3{vars{1:end-1}}),'y-.')
                plot(0,stress(1)-exp1(timecell3{vars{:}}(1)),'rx',0,stress(1)-exp2(timecell3{vars{:}}(1)),'bx')

                %title([astring{alpha_index}, estring{eta_index}, rstring{ramptime_index}, Tstring{T_index} ]); 
                %legend('Data','Single Exp', 'Two Exp', 'Short Time', 'Long Time') 
                %axis([0 1 0 1]);
                xlabel('Time')
                ylabel('Error')
                ylim([-0.1 0.4])
                SaveAsPngEpsAndFig([pwd '/pictures/expfit/relaxation/timestress_diff/' num2str(Tvec(T_index)) '-' num2str(ramptimevec(ramptime_index)) '-' num2str(etavec(eta_index)) '-' num2str(alphavec(alpha_index))]  ,-1, 7, 7/5, 9)
                %pause
                close all
            end
       end
   end
end
%%
%Calculates a three exponential fit. Uses these to plot exponential fits
%and error plots
for T_index = 1
   for ramptime_index = 1:t
       for alpha_index = 1:a
           for  eta_index = 1:g

                vars = {ramptime_index,eta_index,alpha_index,T_index,guess_value};
                figure
                hold on
                exp1 = @(x) fit1(vars{:},1) + fit1(vars{:},2)*exp(-x/fit1(vars{:},3));
                exp2 = @(x) fit2(vars{:},1) + fit2(vars{:},2)*exp(-x/fit2(vars{:},3))+ fit2(vars{:},4)*exp(-x/fit2(vars{:},5));
                exp2_1 = @(x) fit2(vars{:},1)+fit2(vars{:},4)+fit2(vars{:},2)*exp(-x/fit2(vars{:},3));
                exp2_2 = @(x) fit2(vars{:},1) + fit2(vars{:},4)*exp(-x/fit2(vars{:},5));
                exp3 = @(x) three_exp_fit(1) + three_exp_fit(2)*exp(-x/three_exp_fit(3))+ three_exp_fit(4)*exp(-x/three_exp_fit(5))+three_exp_fit(6)*exp(-x/three_exp_fit(7));
 
                time = timecell3{vars{:}}; stress = stresscell3{vars{:}};
                plot(time(1:10:end),stress(1:10:end),'r.',timecell3{vars{:}},exp1(timecell3{vars{:}}),'k-',timecell3{vars{:}},exp2(timecell3{vars{:}}),'b--',timecell3{vars{:}},exp3(timecell3{vars{:}}),'g-.','MarkerSize',4)
                %plot(timecell3{vars{1:end-1}},stresscell3{vars{1:end-1}},'r',timecell3{vars{1:end-1}},exp1(timecell3{vars{1:end-1}}),'k--',timecell3{vars{1:end-1}},exp2(timecell3{vars{1:end-1}}),'b--',timecell3{vars{1:end-1}},exp2_1(timecell3{vars{1:end-1}}),'g-.',timecell3{vars{1:end-1}},exp2_2(timecell3{vars{1:end-1}}),'y-.')
                title([astring{alpha_index}, estring{eta_index}, rstring{ramptime_index}, Tstring{T_index} ]); 
                %legend('Data','Single Exp', 'Two Exp', 'Short Time', 'Long Time') 
                %axis([0 1 0 1]);
                xlabel('Time')
                ylabel('Stress')
                SaveAsPngEpsAndFig(-1,[pwd '/pictures/expfit/creep/timestress_threeexp/' num2str(Tvec(T_index)) '-' num2str(ramptimevec(ramptime_index)) '-' num2str(etavec(eta_index)) '-' num2str(alphavec(alpha_index))]  , 7, 7/5, 9)
                %pause
                figure
                hold on
                plot(timecell3{vars{:}},stress-exp1(timecell3{vars{:}}),'r-',timecell3{vars{:}},stress-exp2(timecell3{vars{:}}),'b-',timecell3{vars{:}},zeros(size(timecell3{vars{:}})),'k--',timecell3{vars{:}},stress-exp3(timecell3{vars{:}}),'g-.')
                plot(0,stress(1)-exp1(timecell3{vars{:}}(1)),'rx',0,stress(1)-exp2(timecell3{vars{:}}(1)),'bx',0,stress(1)-exp3(timecell3{vars{:}}(1)),'gx')

                xlabel('Time')
                ylabel('Error')
                ylim([-0.1 0.35])
                SaveAsPngEpsAndFig([pwd '/pictures/expfit/relaxation/timestress_threeexp_diff/' num2str(Tvec(T_value)) '-' num2str(ramptimevec(ramptime_index)) '-' num2str(etavec(eta_index)) '-' num2str(alphavec(alpha_index))])

            end
       end
   end
end

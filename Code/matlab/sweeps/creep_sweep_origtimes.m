%Varies parameters of stress_2d_ode, and plots time strain graphs in a
%subplot (plus other graphs). Keeps original times without scaling so that
%analysis can be done.
clear all


fbounds = {-1 0 3};
Tbounds = {1 2 2};
etabounds = {-1 0 21};
alphabounds = {-1 0 21};

fstring = arrayfun(@(x) [' Force = 10^{' num2str(x,3) '}'], linspace(fbounds{:}),'UniformOutput',false);
astring = arrayfun(@(x) [' \alpha = 10^{' num2str(x,3) '}'], linspace(alphabounds{:}),'UniformOutput',false);
estring = arrayfun(@(x) [' \eta = 10^{' num2str(x,3) '}'], linspace(etabounds{:}),'UniformOutput',false);
Tstring = arrayfun(@(x) [' T = 10^{' num2str(x,3) '}'], linspace(Tbounds{:}),'UniformOutput',false); Tstring = {' T = 0', Tstring{:}};


forcevec = logspace(fbounds{:});
%Tvec = [0 1 10 100];
Tvec = [0 logspace(Tbounds{:})];
etavec = logspace(etabounds{:});
alphavec = logspace(alphabounds{:});
tend = 200000;
f = length(forcevec);g = length(etavec);a = length(alphavec);
L = f*g*a*length(Tvec);

straincell = cell(1,L);
timecell = cell(1,L);
flagcell = cell(1,L);   
%restoringcell = cell(1,L);
parfor_progress(L);

for i=1:L
    counter = i-1;
    T_index = mod(counter,length(Tvec));
    counter = (counter-T_index)/length(Tvec);
    alpha_index = mod(counter,a);
    counter = (counter-alpha_index)/a;
    eta_index = mod(counter,g);
    counter = (counter-eta_index)/g;
    force_index = counter;
    alpha = alphavec(alpha_index+1);
    force = forcevec(force_index+1);
    eta = etavec(eta_index+1);%converts linear index to alpha,eta,T
    T = Tvec(T_index+1);
    params{i} = {alpha,eta,T,tend,[10,10],force,5};
end
    
parfor i = 1:L%index loops over alpha, then eta, then T
%     counter = i-1;
%     T_index = mod(counter,length(Tvec));
%     counter = (counter-T_index)/length(Tvec);
%     alpha_index = mod(counter,a);
%     counter = (counter-alpha_index)/a;
%     eta_index = mod(counter,g);
%     counter = (counter-eta_index)/g;
%     force_index = counter;
%     alpha = alphavec(alpha_index+1);
%     force = forcevec(force_index+1);
%     eta = etavec(eta_index+1);%converts linear index to alpha,eta,T
%     T = Tvec(T_index+1);
    [Time, Y,~,flag] =stress_2d_ode(params{i}{:});
    N = size(Y,2)/4;
    xvalues = Y(:,1:N);
    strain = (max(xvalues,[],2)-min(xvalues(1,:)))/(max(xvalues(1,:))-min(xvalues(1,:)));
    straincell{i} = strain;
    timecell{i} = Time;
    flagcell{i} = flag;
    parfor_progress;
end
parfor_progress(0);

straincell2 = cell(f,g,a,length(Tvec));
timecell2 = cell(f,g,a,length(Tvec));

fit1 = nan*ones(f,g,a,length(Tvec),3,3);
fit2 = nan*ones(f,g,a,length(Tvec),3,5);
error1 = nan*ones(f,g,a,length(Tvec),3);
error2 = nan*ones(f,g,a,length(Tvec),3);
error_fit_L2 = nan*ones(f,g,a,length(Tvec),3);
error_fit_inf = nan*ones(f,g,a,length(Tvec),3);
%save([pwd '/workspaces/creeperror' num2str(Tvec(1)) '_' num2str(Tvec(end)) '_' num2str(etavec(1)) '-' num2str(etavec(end)) '_' num2str(alphavec(1)) '-' num2str(alphavec(end)) '.mat']);
%%
%load([pwd '/workspaces/creeperror' num2str(Tvec(1)) '_' num2str(Tvec(end)) '_' num2str(etavec(1)) '-' num2str(etavec(end)) '_' num2str(alphavec(1)) '-' num2str(alphavec(end)) '.mat']);
maxstraincell = cell(f,g,a,length(Tvec));
timeendcell = cell(f,g,a,length(Tvec));
straincell3 = cell(f,g,a,length(Tvec));
timecell3 = cell(f,g,a,length(Tvec));
for i = 1:L
    counter = i-1;
    T_index = mod(counter,length(Tvec));
    counter = (counter-T_index)/length(Tvec);
    alpha_index = mod(counter,a);
    counter = (counter-alpha_index)/a;
    eta_index = mod(counter,g);
    counter = (counter-eta_index)/g;
    force_index = counter;
    vars = {force_index+1,eta_index+1,alpha_index+1,T_index+1};
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
        T_index = mod(counter,length(Tvec));
        counter = (counter-T_index)/length(Tvec);
        alpha_index = mod(counter,a);
        counter = (counter-alpha_index)/a;
        eta_index = mod(counter,g);
        counter = (counter-eta_index)/g;
        force_index = counter;
        vars = {force_index+1,eta_index+1,alpha_index+1,T_index+1,guess_index};
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
                for T_index = 1:length(Tvec)
                    vars = {force_index,eta_index,alpha_index,T_index,guess_index};
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
coef_scale_vals = nan*ones(f,g,a,length(Tvec),3);
coef_scale_vals1 = nan*ones(f,g,a,length(Tvec),3);
coef_scale_vals2 = nan*ones(f,g,a,length(Tvec),3);
time_dif_1 =  nan*ones(f,g,a,length(Tvec),3);
time_dif_2 =  nan*ones(f,g,a,length(Tvec),3);
time_dif_3 =  nan*ones(f,g,a,length(Tvec),3);
%Two coefficients of the exponential in the biexponential fit that are the
%same will have a  coef_scale_value of 1, while one which has only 1
%timescale, ie one of the coefficients is zero, will have coef_scale value
%of 0
for guess_index = 1:3
    for force_index = 1:f
        for eta_index = 1:g
            for alpha_index = 1:a;
                for T_index = 1:length(Tvec)
                    vars = {force_index,eta_index,alpha_index,T_index,guess_index};
                    coef1 = abs(fit2(vars{:},2));
                    coef2 = abs(fit2(vars{:},4));
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
clear straincell timecell straincell2 timecell2 %deletes original values
%can take a large amount of space due to long simulation times 
guess_value = 1;
T_value = 1;

viewscale = 11; %Makes subplot display viewscale*viewscale for easier viewing

viewscale = viewscale-1; %This is a poor way to do this, should just use : vector probabaly
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

save([pwd '/workspaces/creeperrortimeorig' num2str(Tvec(1)) '_' num2str(Tvec(end)) '_' num2str(etavec(1)) '-' num2str(etavec(end)) '_' num2str(alphavec(1)) '-' num2str(alphavec(end)) '.mat']);
%% plots strain-time graphs and exponential fits
load([pwd '/workspaces/creeperrortimeorig' num2str(Tvec(1)) '_' num2str(Tvec(end)) '_' num2str(etavec(1)) '-' num2str(etavec(end)) '_' num2str(alphavec(1)) '-' num2str(alphavec(end)) '.mat']);

%%
%Shows example graphs in a set of subplots varying eta and alpha
for force_index = 1:f
    figure
    hold on
    for alpha_index = 1:a_temp 
        for  eta_index = 1:g_temp
            subplot(a_temp,g_temp,eta_index-g_temp*(alpha_index)+a_temp*g_temp)
            vars = {force_index,(eta_index-1)*g_scale+1,(alpha_index-1)*a_scale+1,T_value,guess_value};
            exp1 = @(x) fit1(vars{:},1) + fit1(vars{:},2)*exp(-x/fit1(vars{:},3));
            exp2 = @(x) fit2(vars{:},1) + fit2(vars{:},2)*exp(-x/fit2(vars{:},3))+ fit2(vars{:},4)*exp(-x/fit2(vars{:},5));
            exp2_1 = @(x) fit2(vars{:},1)+fit2(vars{:},4)+fit2(vars{:},2)*exp(-x/fit2(vars{:},3));
            exp2_2 = @(x) fit2(vars{:},1) + fit2(vars{:},4)*exp(-x/fit2(vars{:},5));
            %plot(timecell3{vars{:}},straincell3{vars{:}},'r',timecell3{vars{:}},exp1(timecell3{vars{:}}),'k--',timecell3{vars{:}},exp2(timecell3{vars{:}}),'b--')
            plot(timecell3{vars{1:end-1}},straincell3{vars{1:end-1}},'r',timecell3{vars{1:end-1}},exp1(timecell3{vars{1:end-1}}),'k--',timecell3{vars{1:end-1}},exp2(timecell3{vars{1:end-1}}),'b--',timecell3{vars{1:end-1}},exp2_1(timecell3{vars{1:end-1}}),'g-.',timecell3{vars{1:end-1}},exp2_2(timecell3{vars{1:end-1}}),'y-.')
            title([astring{(alpha_index-1)*a_scale+1}, estring{(eta_index-1)*g_scale+1}, fstring{force_index}, Tstring{T_value} ]); 
            %axis([0 1 0 1]);
            %legend('Data','Single Exp', 'Two Exp', 'Short Time', 'Long Time') 

        end
    end
    %SaveAsPngEpsAndFig(-1,[pwd '/pictures/strainfitting/strainplot-' num2str(T) '-' num2str(timevec(time_index))]  , 7, 7/5, 9)
end

%%
for force_index = 1:f
    for alpha_index = 1:a 
        for  eta_index = 1:g

            %subplot(a_temp,g_temp,eta_index-g_temp*(alpha_index)+a_temp*g_temp)
            %vars = {force_index,(eta_index-1)*g_scale+1,(alpha_index-1)*a_scale+1,T_value,guess_value};
            vars = {force_index,eta_index,alpha_index,T_value,guess_value};

            figure
            hold on
            exp1 = @(x) fit1(vars{:},1) + fit1(vars{:},2)*exp(-x/fit1(vars{:},3));
            exp2 = @(x) fit2(vars{:},1) + fit2(vars{:},2)*exp(-x/fit2(vars{:},3))+ fit2(vars{:},4)*exp(-x/fit2(vars{:},5));
            exp2_1 = @(x) fit2(vars{:},1)+fit2(vars{:},4)+fit2(vars{:},2)*exp(-x/fit2(vars{:},3));
            exp2_2 = @(x) fit2(vars{:},1) + fit2(vars{:},4)*exp(-x/fit2(vars{:},5));
            time = timecell3{vars{1:end-1}}; strain = straincell3{vars{1:end-1}};
            plot(time(1:10:end),strain(1:10:end),'k.',timecell3{vars{:}},exp1(timecell3{vars{:}}),'r-',timecell3{vars{:}},exp2(timecell3{vars{:}}),'b--','MarkerSize',4)

            %plot(timecell3{vars{:}},stresscell3{vars{:}},'r',timecell3{vars{:}},exp1(timecell3{vars{:}}),'k--',timecell3{vars{:}},exp2(timecell3{vars{:}}),'b--')
%             plot(timecell3{vars{1:end-1}},straincell3{vars{1:end-1}},'r',timecell3{vars{1:end-1}},exp1(timecell3{vars{1:end-1}}),'k--',timecell3{vars{1:end-1}},exp2(timecell3{vars{1:end-1}}),'b--',timecell3{vars{1:end-1}},exp2_1(timecell3{vars{1:end-1}}),'g-.',timecell3{vars{1:end-1}},exp2_2(timecell3{vars{1:end-1}}),'y-.')
            %title([astring{(alpha_index-1)*a_scale+1}, estring{(eta_index-1)*g_scale+1}, fstring{force_index}, Tstring{T_value} ]); 

            title([astring{alpha_index}, estring{eta_index}, fstring{force_index}, Tstring{T_value} ]); 

            %legend('Data','Single Exp', 'Two Exp', 'Short Time', 'Long Time') 
           % axis([0 450 0 1]);
            xlabel('Time')
            ylabel('Strain')
            SaveAsPngEpsAndFig([pwd '/pictures/expfit/creep/timestrain/' num2str(Tvec(T_value)) '-' num2str(forcevec(force_index)) '-' num2str(etavec(eta_index)) '-' num2str(alphavec(alpha_index))])
            close all
        end
    end
end
% end
%%
for force_index = 1:f;
    figure
    [X,Y] = meshgrid(etavec,alphavec);
    hold on;
    surf(X,Y,reshape(cell2mat(timeendcell(force_index,:,:,T_value,guess_value)),[g,a])');
    shading interp; 
%     clear alpha
%     alpha(0.5);
    colorbar;
    %contour(X,Y,reshape(cell2mat(timeendcell(force_index,:,:)),[g,a])',contours,'ShowText','on');

    set(gca, 'XScale', 'log', 'YScale', 'log','ZScale','log');
    xlabel('\eta');
    ylabel('\alpha');
    title(['equilibriation times,' fstring{force_index}]);
    SaveAsPngEpsAndFig(-1,[pwd '/pictures/expfit/creep/timeendsurface/' num2str(Tvec(T_value)) '-' num2str(forcevec(force_index))]  , 7, 7/5, 9)
    %SaveAsPngEpsAndFig(-1,[pwd 'asdf']  , 7, 7/5, 9)

end
%%
% for force_index = 1:f
%     for alpha_index = 1:4:a
%         figure
%         hold on
%         for eta_index = 1:g
%             vars = {force_index,eta_index,alpha_index,T_index,guess_value};
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
for force_index = 1:f %creates line plots in order to view differing coefficient values for eta and alpha
    for alpha_index = 1:a_temp
        close all
        figure
        hold on
%         yyaxis right
%         h(1) = plot(etavec,reshape(error1(force_index,:,(alpha_index-1)*a_scale+1,T_value,guess_value),[g 1]),'m-');
%         h(2) = plot(etavec,reshape(error2(force_index,:,(alpha_index-1)*a_scale+1,T_value,guess_value),[g 1]),'g-');
%         ylabel('L2 error')
%         set(gca,'XScale','log','YScale','log')

        for eta_index = 1:g
%             yyaxis left
            vars = {force_index,eta_index,(alpha_index-1)*a_scale+1,T_value,guess_value};
            h(3) = plot(etavec(eta_index),fit2(vars{:},3),'k.','markers',10*coef_scale_vals1(vars{:}));%,'MarkerEdgeColor',(1-coef_scale_vals1(vars{:}))*[1 1 1])
            h(4) = plot(etavec(eta_index),fit2(vars{:},5),'b.','markers',10*coef_scale_vals2(vars{:}));%z`,'MarkerEdgeColor',[1 1 1] + coef_scale_vals2(vars{:})*[-1 -1 0])
            h(5) = plot(etavec(eta_index),fit1(vars{:},3),'r.','markers',10);%*coef_scale_vals2(vars{:}));
%             h(6) = plot(etavec(eta_index),-fit2(vars{:},2),'kx','markers',5);%*coef_scale_vals1(vars{:}))
%             h(7) = plot(etavec(eta_index),-fit2(vars{:},4),'bx','markers',5);%*coef_scale_vals2(vars{:}));
%             h(8) = plot(etavec(eta_index),-fit1(vars{:},2),'rx','markers',5);

            title(['time coefficients,' fstring{force_index} astring{(alpha_index-1)*a_scale+1}])
            set(gca,'XScale','log','YScale','log')

        end
            xlabel('\eta')
            ylabel('Time coefficients')
            SaveAsPngEpsAndFig(-1,[pwd '/pictures/expfit/creep/lineplot/' num2str(Tvec(T_value)) '-' num2str(forcevec(force_index)) '-alpha=' num2str(alphavec((alpha_index-1)*a_scale+1))]  , 7, 7/5, 9)
            %legend(h)
    end
end
%%
for force_index = 1:f %creates line plots in order to view differing coefficient values for eta and alpha
    for eta_index = 1:g_temp
        close all
        figure
        hold on
%         yyaxis right
%         h(1) = plot(alphavec,reshape(error1(force_index,(eta_index-1)*g_scale+1,:,T_value,guess_value),[a 1]),'m-');
%         h(2) = plot(alphavec,reshape(error2(force_index,(eta_index-1)*g_scale+1,:,T_value,guess_value),[a 1]),'g-');
%         ylabel('L2 error')
%         set(gca,'XScale','log','YScale','log')

        for alpha_index = 1:a
%             yyaxis left
            vars = {force_index,(eta_index-1)*g_scale+1,alpha_index,T_value,guess_value};
            h(3) = plot(alphavec(alpha_index),fit2(vars{:},3),'k.','markers',10*coef_scale_vals1(vars{:}));%,'MarkerEdgeColor',(1-coef_scale_vals1(vars{:}))*[1 1 1])
            h(4) = plot(alphavec(alpha_index),fit2(vars{:},5),'b.','markers',10*coef_scale_vals2(vars{:}));%z`,'MarkerEdgeColor',[1 1 1] + coef_scale_vals2(vars{:})*[-1 -1 0])
            h(5) = plot(alphavec(alpha_index),fit1(vars{:},3),'r.','markers',10);%*coef_scale_vals2(vars{:}));
%             h(6) = plot(alphavec(alpha_index),-fit2(vars{:},2),'kx','markers',5);%*coef_scale_vals1(vars{:}))
%             h(7) = plot(alphavec(alpha_index),-fit2(vars{:},4),'bx','markers',5);%*coef_scale_vals2(vars{:}));
%             h(8) = plot(alphavec(alpha_index),-fit1(vars{:},2),'rx','markers',5);

            title(['time coefficients,' fstring{force_index} estring{(eta_index-1)*g_scale+1}])
            set(gca,'XScale','log','YScale','log')

        end
            xlabel('\alpha')
            ylabel('Time coefficients')
            SaveAsPngEpsAndFig(-1,[pwd '/pictures/expfit/creep/lineplot/' num2str(Tvec(T_value)) '-' num2str(forcevec(force_index)) '-eta=' num2str(etavec((eta_index-1)*g_scale+1))]  , 7, 7/5, 9)
            %legend(h)
    end
end
%%
[X,Y] = meshgrid(etavec,alphavec);
for force_index = 1:f
    figure
    surf(X,Y,reshape(fit2(force_index,:,:,T_value,guess_value,3),[g,a])')
    xlabel('\eta');
    ylabel('\alpha');
    title(['Timescale 1,' fstring{force_index}])
    set(gca, 'XScale', 'log', 'YScale', 'log');
    colorbar;
end

for force_index = 1:f
    figure
    surf(X,Y,reshape(fit2(force_index,:,:,T_value,guess_value,5),[g,a])')
    xlabel('\eta');
    ylabel('\alpha');
    title(['Timescale 2,' fstring{force_index}])
    set(gca, 'XScale', 'log', 'YScale', 'log');
    colorbar;
end

%%
[X,Y] = meshgrid(etavec,alphavec);

% for force_index = 1:f
%     figure
%     surf(X,Y,reshape(time_dif_1(force_index,:,:,T_value,guess_value),[g,a])')
%     xlabel('eta');
%     ylabel('alpha');
%     title(['Timedif1, force = ' num2str(forcevec(force_index))])
%     set(gca, 'XScale', 'log', 'YScale', 'log');
%     colorbar;
% end

for force_index = 1:f
    figure
    surf(X,Y,reshape(time_dif_2(force_index,:,:,T_value,guess_value),[g,a])')
    xlabel('\eta');
    ylabel('\alpha');
    title(['CoefRatio,' fstring{force_index}])
    set(gca, 'XScale', 'log', 'YScale', 'log');
    colorbar;
end

% for force_index = 1:f
%     figure
%     surf(X,Y,reshape(time_dif_3(force_index,:,:,T_value,guess_value),[g,a])')
%     xlabel('eta');
%     ylabel('alpha');
%     title(['Timedif3, force = ' num2str(forcevec(force_index))])
%     set(gca, 'XScale', 'log', 'YScale', 'log');
%     colorbar;
% end

%%
[X,Y] = meshgrid(etavec,alphavec);

for force_index = 1:f
    figure
    surf(X,Y,reshape(error1(force_index,:,:,T_value,guess_value),[g,a])')
    
    xlabel('\eta');
    ylabel('\alpha');
    title(['one exponential L2 error,' fstring{force_index}])
    set(gca, 'XScale', 'log', 'YScale', 'log');
end

%%

[X,Y] = meshgrid(etavec,alphavec);
time_dif_2(:,:,a,:,:) = nan*time_dif_2(:,:,a,:,:);

contours = logspace(0,5,61);
for force_index = 1:f %view the difference between coefficients for one and two exp fits as a surface plot
    for T_index = 1:length(Tvec)
        figure
        hold on
        surf(X,Y,reshape(time_dif_2(force_index,:,:,T_index,guess_value),[g,a])')
        clearvars alpha
        %alpha(0.5)
        %contour(X,Y,reshape(time_dif_2(force_index,:,:,T_index,guess_value),[g,a])',contours,'ShowText','on')
        xlabel('\eta');
        ylabel('\alpha');
        title(['CoefRatio,' fstring{force_index} Tstring{T_index}])
        set(gca, 'XScale', 'log', 'YScale', 'log');
        colorbar;
        caxis([1,3])
        SaveAsPngEpsAndFig(-1,[pwd '/pictures/expfit/creep/timedifsurface/' num2str(Tvec(T_index)) '-' num2str(forcevec(force_index))]  , 7, 7/5, 9)
    end
end

%%
time_dif_2(:,:,a,:,:) = nan*time_dif_2(:,:,a,:,:);
stylevec = {'-','--','-.',':'};
colourvec = {[1 0 0],[0 0 1],[0 0 0],[0 1 0],[1 1 0],[1 0 1],[0 1 1],[1 0 0],[0 0 1],[0 0 0],[0 1 0]};%need to stop reuse of colours, temp measure
error_threshold = 1.5;
figure
hold on;
[X,Y] = meshgrid(etavec,alphavec);

two_exp_status = time_dif_2 > error_threshold;

for T_index = 1:length(Tvec) %view compound contour plot - generally too messy to use
    for force_index = 1:f
        [~,h((T_index-1)*f+force_index)] = contour(X,Y,reshape(time_dif_2(force_index,:,:,T_index,guess_value),[g,a])',[error_threshold,error_threshold],[stylevec{T_index}],'ShowText','off','LineColor',colourvec{force_index});
        set(gca, 'XScale', 'log', 'YScale', 'log');
        legendcell{((T_index-1)*f+force_index)} =['T = ', num2str(Tvec(T_index)), ',force = ' , num2str(forcevec(force_index))];
        xlabel('\eta');
        ylabel('\alpha');
        title(['overall relaxation contour, threshold = ' num2str(error_threshold)])

    end
end

legendflex(h,legendcell)

%%
clear legendcell h
colourvec = {'r','b','g'};
sizevec = {4,5,6};
for T_index = 1:1%3:length(Tvec) %v
    figure
    hold on
    clear h legendcell 
    %view regions of two exponential bejaviour for varying force
    for force_index = f:-1:1
        [~,h(force_index)] = contour(X,Y,reshape(time_dif_2(force_index,:,:,T_index,guess_value),[g,a])',[error_threshold,error_threshold],[stylevec{1}],'ShowText','off','LineColor',colourvec{force_index});
        set(gca, 'XScale', 'log', 'YScale', 'log');
        %legendcell{(force_index)} =['force = ' , num2str(forcevec(force_index))];
        xlabel('\eta');
        ylabel('\alpha');
        title([%'overall relaxation contour, 
            'Threshold = ' num2str(error_threshold), Tstring{T_index}])
        Z = reshape(two_exp_status(force_index,:,:,T_index,guess_value),[g,a])';
        plot(Z.*X,Z.*Y,['.'],'markers',sizevec{force_index},'Color',colourvec{force_index})
    end
    legendflex(h,fstring)
    axis([0.1 1 0.1 1])

    SaveAsPngEpsAndFig([pwd '/pictures/expfit/creep/biexpcontour/' num2str(Tvec(T_index))])


end


%%
clear legendcell h
for force_index = 2:2%1:f
    figure
    hold on
    clear h legendcell
        %view regions of two exponential bejaviour for varying T_mem
        %averaging values
    for T_index = length(Tvec):-1:1
        [~,h(T_index)] = contour(X,Y,reshape(time_dif_2(force_index,:,:,T_index,guess_value),[g,a])',[error_threshold,error_threshold],[stylevec{1}],'ShowText','off','LineColor',colourvec{T_index});
        set(gca, 'XScale', 'log', 'YScale', 'log');
        %legendcell{(T_index)} =['T = ' , num2str(Tvec(T_index))];
        xlabel('\eta');
        ylabel('\alpha');
        title(['Threshold = ' num2str(error_threshold), fstring{force_index}])
        Z = reshape(two_exp_status(force_index,:,:,T_index,guess_value),[g,a])';
        plot(Z.*X,Z.*Y,['.'],'markers',sizevec{T_index},'Color',colourvec{T_index} )
        %surf(Z.*X,Z.*Y,+Z)
        %colormap(colourvec{T_index})
        %alpha(0.3)
        %shading interp
    end
%     h2 = [h(1) h(4) h(6)];
%     Tstring2 ={Tstring{1} Tstring{4} Tstring{6}};
    legendflex(h,Tstring)
    axis([0.1 1 0.1 1])

    %SaveAsPngEpsAndFig([pwd '/pictures/expfit/creep/biexpcontour2/' num2str(forcevec(force_index))])

end


%%
% vars = {force_index,eta_index,alpha_index,T_value,guess_value};
% %[~,~,~,~,~,~,three_exp_fit] = CalculateExponentialFits(timecell3{vars{1:end-1}}',straincell3{vars{1:end-1}}',1,guess_index);
% 
% 
% % figure
% % hold on
% % %subplot(a_temp,g_temp,eta_index-g_temp*(alpha_index)+a_temp*g_temp)
% exp1 = @(x) fit1(vars{:},1) + fit1(vars{:},2)*exp(-x/fit1(vars{:},3));
% exp2 = @(x) fit2(vars{:},1) + fit2(vars{:},2)*exp(-x/fit2(vars{:},3))+ fit2(vars{:},4)*exp(-x/fit2(vars{:},5));
% exp2_1 = @(x) fit2(vars{:},1)+fit2(vars{:},4)+fit2(vars{:},2)*exp(-x/fit2(vars{:},3));
% exp2_2 = @(x) fit2(vars{:},1) + fit2(vars{:},4)*exp(-x/fit2(vars{:},5));
% %exp3 = @(x) three_exp_fit(1) + three_exp_fit(2)*exp(-x/three_exp_fit(3))+ three_exp_fit(4)*exp(-x/three_exp_fit(5))+three_exp_fit(6)*exp(-x/three_exp_fit(7));
% 
% % time = timecell3{vars{:}}; strain = straincell3{vars{:}};
% % %plot(time(1:10:end),strain(1:10:end),'r.',timecell3{vars{:}},exp1(timecell3{vars{:}}),'k-',timecell3{vars{:}},exp2(timecell3{vars{:}}),'b--',timecell3{vars{:}},exp3(timecell3{vars{:}}),'g-.','MarkerSize',4)
% % %plot(timecell3{vars{1:end-1}},straincell3{vars{1:end-1}},'r',timecell3{vars{1:end-1}},exp1(timecell3{vars{1:end-1}}),'k--',timecell3{vars{1:end-1}},exp2(timecell3{vars{1:end-1}}),'b--',timecell3{vars{1:end-1}},exp2_1(timecell3{vars{1:end-1}}),'g-.',timecell3{vars{1:end-1}},exp2_2(timecell3{vars{1:end-1}}),'y-.')
% % %title([astring{alpha_index}, estring{eta_index}, fstring{force_index}, Tstring{T_value} ]); 
% % %legend('Data','Single Exp', 'Two Exp', 'Short Time', 'Long Time') 
% % %axis([0 1 0 1]);
% % xlabel('Time')
% % ylabel('Strain')
% %SaveAsPngEpsAndFig(-1,[pwd '/pictures/expfit/creep/timestrain_threeexp/' num2str(Tvec(T_value)) '-' num2str(forcevec(force_index)) '-' num2str(etavec(eta_index)) '-' num2str(alphavec(alpha_index))]  , 7, 7/5, 9)
% %pause
% 
% figure
% hold on
% %subplot(a_temp,g_temp,eta_index-g_temp*(alpha_index)+a_temp*g_temp)
% exp1 = @(x) fit1(vars{:},1) + fit1(vars{:},2)*exp(-x/fit1(vars{:},3));
% exp2 = @(x) fit2(vars{:},1) + fit2(vars{:},2)*exp(-x/fit2(vars{:},3))+ fit2(vars{:},4)*exp(-x/fit2(vars{:},5));
% exp2_1 = @(x) fit2(vars{:},1)+fit2(vars{:},4)+fit2(vars{:},2)*exp(-x/fit2(vars{:},3));
% exp2_2 = @(x) fit2(vars{:},1) + fit2(vars{:},4)*exp(-x/fit2(vars{:},5));
% %exp3 = @(x) three_exp_fit(1) + three_exp_fit(2)*exp(-x/three_exp_fit(3))+ three_exp_fit(4)*exp(-x/three_exp_fit(5))+three_exp_fit(6)*exp(-x/three_exp_fit(7));
% time = timecell3{vars{:}}; strain = straincell3{vars{:}};
% %plot(timecell3{vars{:}},strain-exp1(timecell3{vars{:}}),'k-',timecell3{vars{:}},strain-exp2(timecell3{vars{:}}),'b--',timecell3{vars{:}},strain-exp3(timecell3{vars{:}}),'g-.',timecell3{vars{:}},zeros(size(timecell3{vars{:}})),'r--','MarkerSize',4)
% %plot(timecell3{vars{:}},strain-exp3(timecell3{vars{:}}),'g-',timecell3{vars{:}},strain-exp2(timecell3{vars{:}}),'b--',timecell3{vars{:}},zeros(size(timecell3{vars{:}})),'r--','MarkerSize',4)
% plot(timecell3{vars{:}},strain-exp1(timecell3{vars{:}}),'k-',timecell3{vars{:}},strain-exp2(timecell3{vars{:}}),'b',timecell3{vars{:}},zeros(size(timecell3{vars{:}})),'r--')%,timecell3{vars{:}},strain-exp3(timecell3{vars{:}}),'g-.','MarkerSize',4)  
% %plot(timecell3{vars{1:end-1}},straincell3{vars{1:end-1}},'r',timecell3{vars{1:end-1}},exp1(timecell3{vars{1:end-1}}),'k--',timecell3{vars{1:end-1}},exp2(timecell3{vars{1:end-1}}),'b--',timecell3{vars{1:end-1}},exp2_1(timecell3{vars{1:end-1}}),'g-.',timecell3{vars{1:end-1}},exp2_2(timecell3{vars{1:end-1}}),'y-.')
% 'asdf'
% 
% 
% %plot(timecell3{vars{1:end-1}},straincell3{vars{1:end-1}},'r',timecell3{vars{1:end-1}},exp1(timecell3{vars{1:end-1}}),'k--',timecell3{vars{1:end-1}},exp2(timecell3{vars{1:end-1}}),'b--',timecell3{vars{1:end-1}},exp2_1(timecell3{vars{1:end-1}}),'g-.',timecell3{vars{1:end-1}},exp2_2(timecell3{vars{1:end-1}}),'y-.')
% %title([astring{alpha_index}, estring{eta_index}, fstring{force_index}, Tstring{T_value} ]); 
% %legend('Data','Single Exp', 'Two Exp', 'Short Time', 'Long Time') 
% %axis([0 1 0 1]);
% xlabel('Time')
% ylabel('Error')
% ylim([-0.2,0.05])
% 
% %SaveAsPngEpsAndFig(-1,['/Users/reubenv/Desktop/expfit/creep/timestrain_diff/' num2str(Tvec(T_value)) '-' num2str(forcevec(force_index)) '-' num2str(etavec(eta_index)) '-' num2str(alphavec(alpha_index))]  , 7, 7/5, 9)
% %pause
% figure
% plot(timecell3{vars{1:end-1}},straincell3{vars{1:end-1}},'r',timecell3{vars{1:end-1}},exp1(timecell3{vars{1:end-1}}),'k--',timecell3{vars{1:end-1}},exp2(timecell3{vars{1:end-1}}),'b--',timecell3{vars{1:end-1}},exp2_1(timecell3{vars{1:end-1}}),'g-.',timecell3{vars{1:end-1}},exp2_2(timecell3{vars{1:end-1}}),'y-.')
% legend('Data','Single Exp', 'Two Exp', 'Short Time', 'Long Time') 
% xlabel('Time')
% ylabel('Strain')
%% Calculates error plots for one and two exp fits
for force_index = 1:f
    for alpha_index = 1:a
        for  eta_index = 1:g
            for T_index = 1:T
                vars = {force_index,eta_index,alpha_index,T_value,guess_value};
                figure
                hold on
                %subplot(a_temp,g_temp,eta_index-g_temp*(alpha_index)+a_temp*g_temp)
                exp1 = @(x) fit1(vars{:},1) + fit1(vars{:},2)*exp(-x/fit1(vars{:},3));
                exp2 = @(x) fit2(vars{:},1) + fit2(vars{:},2)*exp(-x/fit2(vars{:},3))+ fit2(vars{:},4)*exp(-x/fit2(vars{:},5));
                exp2_1 = @(x) fit2(vars{:},1)+fit2(vars{:},4)+fit2(vars{:},2)*exp(-x/fit2(vars{:},3));
                exp2_2 = @(x) fit2(vars{:},1) + fit2(vars{:},4)*exp(-x/fit2(vars{:},5));
                exp3 = @(x) three_exp_fit(1) + three_exp_fit(2)*exp(-x/three_exp_fit(3))+ three_exp_fit(4)*exp(-x/three_exp_fit(5))+three_exp_fit(6)*exp(-x/three_exp_fit(7));

                time = timecell3{vars{:}}; strain = straincell3{vars{:}};
                hold on
                plot(timecell3{vars{:}},strain-exp1(timecell3{vars{:}}),'r-',timecell3{vars{:}},strain-exp2(timecell3{vars{:}}),'b-',timecell3{vars{:}},zeros(size(timecell3{vars{:}})),'k--')%,timecell3{vars{:}},strain-exp3(timecell3{vars{:}}),'g-.','MarkerSize',4)            %plot(timecell3{vars{1:end-1}},straincell3{vars{1:end-1}},'r',timecell3{vars{1:end-1}},exp1(timecell3{vars{1:end-1}}),'k--',timecell3{vars{1:end-1}},exp2(timecell3{vars{1:end-1}}),'b--',timecell3{vars{1:end-1}},exp2_1(timecell3{vars{1:end-1}}),'g-.',timecell3{vars{1:end-1}},exp2_2(timecell3{vars{1:end-1}}),'y-.')
                plot(0,strain(1)-exp1(timecell3{vars{:}}(1)),'rx',0,strain(1)-exp2(timecell3{vars{:}}(1)),'bx')
                xlabel('Time')
                ylabel('Error')
                ylim([-0.2,0.05])
                SaveAsPngEpsAndFig([pwd '/pictures/expfit/creep/timestrain_diff/' num2str(Tvec(T_value)) '-' num2str(forcevec(force_index)) '-' num2str(etavec(eta_index)) '-' num2str(alphavec(alpha_index))])
                %pause
                close all
            end
        end
    end
end
%%
%Calculates a three exponential fit. Uses these to plot exponential fits
%and error plots
for force_index = 1:f
    for alpha_index = 1:a
        for  eta_index = 1:g
            for T_index = 1:T
         %   figure
         %   hold on
            %subplot(a_temp,g_temp,eta_index-g_temp*(alpha_index)+a_temp*g_temp)
            vars = {force_index,eta_index,alpha_index,T_value,guess_value};
               figure
               hold on
                subplot(a_temp,g_temp,eta_index-g_temp*(alpha_index)+a_temp*g_temp)
                exp1 = @(x) fit1(vars{:},1) + fit1(vars{:},2)*exp(-x/fit1(vars{:},3));
                exp2 = @(x) fit2(vars{:},1) + fit2(vars{:},2)*exp(-x/fit2(vars{:},3))+ fit2(vars{:},4)*exp(-x/fit2(vars{:},5));
                exp2_1 = @(x) fit2(vars{:},1)+fit2(vars{:},4)+fit2(vars{:},2)*exp(-x/fit2(vars{:},3));
                exp2_2 = @(x) fit2(vars{:},1) + fit2(vars{:},4)*exp(-x/fit2(vars{:},5));
               exp3 = @(x) three_exp_fit(1) + three_exp_fit(2)*exp(-x/three_exp_fit(3))+ three_exp_fit(4)*exp(-x/three_exp_fit(5))+three_exp_fit(6)*exp(-x/three_exp_fit(7));

                time = timecell3{vars{:}}; strain = straincell3{vars{:}};
                plot(time(1:10:end),strain(1:10:end),'r.',timecell3{vars{:}},exp1(timecell3{vars{:}}),'k-',timecell3{vars{:}},exp2(timecell3{vars{:}}),'b--',timecell3{vars{:}},exp3(timecell3{vars{:}}),'g-.','MarkerSize',4)
                %plot(timecell3{vars{1:end-1}},straincell3{vars{1:end-1}},'r',timecell3{vars{1:end-1}},exp1(timecell3{vars{1:end-1}}),'k--',timecell3{vars{1:end-1}},exp2(timecell3{vars{1:end-1}}),'b--',timecell3{vars{1:end-1}},exp2_1(timecell3{vars{1:end-1}}),'g-.',timecell3{vars{1:end-1}},exp2_2(timecell3{vars{1:end-1}}),'y-.')
                title([astring{alpha_index}, estring{eta_index}, fstring{force_index}, Tstring{T_value} ]); 
                %legend('Data','Single Exp', 'Two Exp', 'Short Time', 'Long Time') 
                %axis([0 1 0 1]);
                xlabel('Time')
                ylabel('Strain')
                SaveAsPngEpsAndFig(-1,[pwd '/pictures/expfit/creep/timestrain_threeexp/' num2str(Tvec(T_value)) '-' num2str(forcevec(force_index)) '-' num2str(etavec(eta_index)) '-' num2str(alphavec(alpha_index))]  , 7, 7/5, 9)
                %pause
                figure
                hold on
                plot(timecell3{vars{:}},strain-exp1(timecell3{vars{:}}),'k-',timecell3{vars{:}},strain-exp2(timecell3{vars{:}}),'b',timecell3{vars{:}},zeros(size(timecell3{vars{:}})),'r--',timecell3{vars{:}},strain-exp3(timecell3{vars{:}}),'g-.','MarkerSize',4)  
                title([astring{alpha_index}, estring{eta_index}, fstring{force_index}, Tstring{T_value} ]); 

                xlabel('Time')
                ylabel('Error')
                SaveAsPngEpsAndFig(-1,[pwd '/pictures/expfit/creep/timestrain_threeexp_diff/' num2str(Tvec(T_value)) '-' num2str(forcevec(force_index)) '-' num2str(etavec(eta_index)) '-' num2str(alphavec(alpha_index))]  , 7, 7/5, 9)

            end
        end
    end
end


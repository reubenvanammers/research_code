%Varies parameters of stress_2d_ode, and plots time strain graphs in a
%subplot.
clear all

forcevec = logspace(-1.5,-1 ,3);
%Tvec = [0 20 100];
T = 0;
etavec = logspace(-2,0,17);
alphavec = logspace(-2,0,17);
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

fit1 = nan*ones(f,g,a,3);
fit2 = nan*ones(f,g,a,5);
error1 = nan*ones(f,g,a);
error2 = nan*ones(f,g,a);
error_fit_L2 = nan*ones(f,g,a);
error_fit_inf = nan*ones(f,g,a);
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
    time_scaled = timecell2{vars{:}}/timeendcell{vars{:}};
    
    timecell3{vars{:}} = linspace(x_coord_scale(1),x_coord_scale(2),1001)';
    straincell3{vars{:}} = interp1(time_scaled,strain_scaled,timecell3{vars{:}});
%     timecell3{vars{:}} = linspace(0,timeendcell{vars{:}},1001)';
%     straincell3{vars{:}} = interp1(timecell2{vars{:}},straincell2{vars{:}},timecell3{vars{:}});

end
%%
unable_to_fit = [];
for i = 1:L
    counter = i-1;
    alpha_index = mod(counter,a);
    counter = (counter-alpha_index)/a;
    eta_index = mod(counter,g);
    counter = (counter-eta_index)/g;
    force_index = counter;
    vars = {force_index+1,eta_index+1,alpha_index+1};
%     straincell2{vars{:}} = straincell{i};
%     timecell2{vars{:}} = timecell{i};
    if ~flagcell{i}
        try
            [fit1(vars{:},:),error1(vars{:}),fit2(vars{:},:),error2(vars{:}),error_fit_L2(vars{:}),error_fit_inf(vars{:})] = CalculateExponentialFits(timecell3{vars{:}}',straincell3{vars{:}}');
        catch
            unable_to_fit = [unable_to_fit; vars];
        end
    else
        %unable_to_fit = [unable_to_fit; vars];
    end
end
clear straincell timecell straincell2 timecell2
save([pwd '/workspaces/creeperror' num2str(T) '_' num2str(etavec(1)) '-' num2str(etavec(end)) '_' num2str(alphavec(1)) '-' num2str(alphavec(end)) '.mat']);
%% plots strain-time graphs and exponential fits
load([pwd '/workspaces/creeperror' num2str(T) '_' num2str(etavec(1)) '-' num2str(etavec(end)) '_' num2str(alphavec(1)) '-' num2str(alphavec(end)) '.mat']);
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


for force_index = 1:f
    figure
    hold on
    for alpha_index = 1:a_temp 
        for  eta_index = 1:g_temp
            subplot(a_temp,g_temp,eta_index-g_temp*(alpha_index)+a_temp*g_temp)
            vars = {force_index,(eta_index-1)*g_scale+1,(alpha_index-1)*a_scale+1};
            exp1 = @(x) fit1(vars{:},1) + fit1(vars{:},2)*exp(-x/fit1(vars{:},3));
            exp2 = @(x) fit2(vars{:},1) + fit2(vars{:},2)*exp(-x/fit2(vars{:},3))+ fit2(vars{:},4)*exp(-x/fit2(vars{:},5));
            exp2_1 = @(x) fit2(vars{:},1)+fit2(vars{:},4)+fit2(vars{:},2)*exp(-x/fit2(vars{:},3));
            exp2_2 = @(x) fit2(vars{:},1) + fit2(vars{:},4)*exp(-x/fit2(vars{:},5));

%             plot(timecell3{vars{:}},straincell3{vars{:}},'r',timecell3{vars{:}},exp1(timecell3{vars{:}}),'k--',timecell3{vars{:}},exp2(timecell3{vars{:}}),'b--')
            plot(timecell3{vars{:}},straincell3{vars{:}},'r',timecell3{vars{:}},exp1(timecell3{vars{:}}),'k--',timecell3{vars{:}},exp2(timecell3{vars{:}}),'b--',timecell3{vars{:}},exp2_1(timecell3{vars{:}}),'g--',timecell3{vars{:}},exp2_2(timecell3{vars{:}}),'y--')

            title(['alpha = ', num2str(alphavec_temp(alpha_index)), ' eta = ', num2str(etavec_temp(eta_index)), ' force = ' num2str(forcevec(force_index))]); 
        end
    end
%    SaveAsPngEpsAndFig(-1,[pwd '/pictures/creepfitting/strainplot-' num2str(T) '-' strrep(num2str(forcevec(force_index)),'.','')]  , 7, 7/5, 9)
end

%%
error_threshold = 1e-5;
colourvec = ['y','m','c','r','b','g','k'];
figure
hold on
[X,Y] = meshgrid(etavec,alphavec);
contours = logspace(-10,10,21);
plot(X,Y,'ko');
h = [];

for force_index = 1:f
    [~,h(force_index)] = contour(X,Y,reshape(error_fit_L2(force_index,:,:),[g,a])',[error_threshold,error_threshold],colourvec(force_index),'ShowText','off');
    set(gca, 'XScale', 'log', 'YScale', 'log');
end
xlabel('eta');
ylabel('alpha');
title(['regions above/right of contours are good fits/have 1 timescale, have error less than ' num2str(error_threshold)])
legendflex(h,arrayfun(@num2str,forcevec,'UniformOutput',false));
SaveAsPngEpsAndFig(-1,[pwd '/pictures/creepfitting/fitcontour-' num2str(T)]  , 7, 7/5, 9)


%%
contours = logspace(-10,10,21);
for force_index = 1:f;
    figure
    [X,Y] = meshgrid(etavec,alphavec);
    hold on;
    plot(X,Y,'k.','MarkerSize',12)
    surf(X,Y,reshape(error_fit_L2(force_index,:,:),[g,a])');
     shading interp;
     alpha(0.5);
     colorbar;
%     contour(X,Y,reshape(error_fit_inf(force_index,:,:),[g,a])',contours,'ShowText','on');

    set(gca, 'XScale', 'log', 'YScale', 'log');
    xlabel('eta');
    ylabel('alpha');
    title(['Error contour - force = ' num2str(forcevec(force_index))])
%    SaveAsPngEpsAndFig(-1,[pwd '/pictures/creepfitting/forcefitcontour-' num2str(T) '-' strrep(num2str(forcevec(force_index)),'.','')]  , 7, 7/5, 9)

end

%%
for force_index = 1:1;
    figure
    [X,Y] = meshgrid(etavec,alphavec);
    hold on;
    surf(X,Y,reshape(cell2mat(timeendcell(force_index,:,:)),[g,a])');
    shading interp;
    alpha(0.5);
    colorbar;
    %contour(X,Y,reshape(cell2mat(timeendcell(force_index,:,:)),[g,a])',contours,'ShowText','on');

    set(gca, 'XScale', 'log', 'YScale', 'log','ZScale','log');
    xlabel('eta');
    ylabel('alpha');
    title(['equilibriation times, force = ', num2str(forcevec(force_index))]);
    SaveAsPngEpsAndFig(-1,[pwd '/pictures/creepfitting/equilibriationtimes-' num2str(T) '-' strrep(num2str(forcevec(force_index)),'.','')]  , 7, 7/5, 9)
    %SaveAsPngEpsAndFig(-1,[pwd 'asdf']  , 7, 7/5, 9)

end
%%
for force_index = 1:1;
    figure
    [X,Y] = meshgrid(etavec,alphavec);
    hold on;
    surf(X,Y,reshape(cell2mat(maxstraincell(force_index,:,:,T_value,guess_value)),[g,a])');
    shading interp;
    colorbar;
    %contour(X,Y,reshape(cell2mat(maxstraincell(force_index,:,:)),[g,a])',contours,'ShowText','on');

    set(gca, 'XScale', 'log', 'YScale', 'log','ZScale','log');
    xlabel('\eta');
    ylabel('\alpha');
    %title('Max Strain')
    %view([90 0])

    %SaveAsPngEpsAndFig(-1,[pwd '/pictures/creepfitting/equilibriationlenghs-' num2str(T) '-' strrep(num2str(forcevec(force_index)),'.','')]  , 7, 7/5, 9)
end

%%
for force_index = 1:1
    for alpha_index = 1:1
        figure
        error_vals = error_fit_inf(force_index,:,alpha_index);
        semilogx(etavec,error_vals);
        xlabel('eta')
        ylabel('Infinity norm error')
        title(['force = ' num2str(forcevec(force_index)) ', alpha = ' num2str(alphavec(alpha_index))])
    end
end

%%
biexponential_status = nan*zeros(f,g,a);
scale_diff_vals = nan*zeros(f,g,a);
timescale_diff_vals =nan*zeros(f,g,a);
timescale_range = [0.2 5];
scale_range = [0.2 5];
[X,Y] = meshgrid(etavec,alphavec);

for force_index = 1:f
    for eta_index = 1:g
        for alpha_index = 1:a;
            vars = {force_index,eta_index,alpha_index};
            timescale_diff = fit2(vars{:},5)/fit2(vars{:},3);
            timescale_diff_vals(vars{:}) = timescale_diff;
            scale_diff = fit2(vars{:},4)/fit2(vars{:},2);
            scale_diff_vals(vars{:}) = scale_diff;
            if timescale_diff < timescale_range(1) || timescale_diff > timescale_range(2)
                if scale_diff > scale_range(1) && scale_diff < scale_range(2)
                    biexponential_status(vars{:}) =1;
                else-
                    biexponential_status(vars{:}) = 0;
                end
            end
        end
    end
end


for force_index = 1:1
    figure
    surf(X,Y,reshape(biexponential_status(force_index,:,:),[g,a])')
    xlabel('eta');
    ylabel('alpha');
    title(['biexponential status, force = ' num2str(forcevec(force_index))])
    set(gca, 'XScale', 'log', 'YScale', 'log');
    colorbar;
end


for force_index = 1:1
    figure
    surf(X,Y,reshape(scale_diff_vals(force_index,:,:),[g,a])')
    xlabel('eta');
    ylabel('alpha');
    title(['Coefficient scaling difference, force = ' num2str(forcevec(force_index))])
    set(gca, 'XScale', 'log', 'YScale', 'log');
    colorbar;
end

for force_index = 1:1
    figure
    surf(X,Y,reshape(timescale_diff_vals(force_index,:,:),[g,a])')
    xlabel('eta');
    ylabel('alpha');
    title(['exponential timescale difference, force = ' num2str(forcevec(force_index))])
    set(gca, 'XScale', 'log', 'YScale', 'log');
    colorbar;
end


%%
one_exp_error = zeros(f,g,a);
two_exp_error = zeros(f,g,a);
[X,Y] = meshgrid(etavec,alphavec);
for force_index = 1:f %calculates the error at x = 0 (which will also be the maximum error) for each force, eta, and alpha for one and two exponential fits
    for alpha_index = 1:a 
        for  eta_index = 1:g
            vars = {force_index,eta_index,alpha_index};
            exp1 = @(x) fit1(vars{:},1) + fit1(vars{:},2)*exp(-x/fit1(vars{:},3));
            exp2 = @(x) fit2(vars{:},1) + fit2(vars{:},2)*exp(-x/fit2(vars{:},3))+ fit2(vars{:},4)*exp(-x/fit2(vars{:},5));
            one_exp_error(force_index,eta_index,alpha_index) = abs(exp1(0));
            two_exp_error(force_index,eta_index,alpha_index) = abs(exp2(0));
        end
    end
end


for force_index = 1:f
    figure
    surf(X,Y,reshape(one_exp_error(force_index,:,:),[g,a])')
    
    xlabel('eta');
    ylabel('alpha');
    title(['one exponential maximum error, force = ' num2str(forcevec(force_index))])
    set(gca, 'XScale', 'log', 'YScale', 'log');
end


for force_index = 1:f
    figure
    surf(X,Y,reshape(two_exp_error(force_index,:,:),[g,a])')
    
    xlabel('eta');
    ylabel('alpha');
    title(['two exponential maximum error, force = ' num2str(forcevec(force_index))])
    set(gca, 'XScale', 'log', 'YScale', 'log');
end

%%

for force_index = 1:f
    figure
    surf(X,Y,reshape(error1(force_index,:,:),[g,a])')
    
    xlabel('eta');
    ylabel('alpha');
    title(['one exponential L2 error, force = ' num2str(forcevec(force_index))])
    set(gca, 'XScale', 'log', 'YScale', 'log');
end


for force_index = 1:f
    figure
    surf(X,Y,reshape(error2(force_index,:,:),[g,a])')
    
    xlabel('eta');
    ylabel('alpha');
    title(['two exponential L2 error, force = ' num2str(forcevec(force_index))])
    set(gca, 'XScale', 'log', 'YScale', 'log');
end
%%one_exp_error = zeros(f,g,a);
two_exp_error = zeros(f,g,a);
[X,Y] = meshgrid(etavec,alphavec);
for force_index = 1:f %calculates the error at x = 0 (which will also be the maximum error) for each force, eta, and alpha for one and two exponential fits
    for alpha_index = 1:a 
        for  eta_index = 1:g
            vars = {force_index,eta_index,alpha_index};
            exp1 = @(x) fit1(vars{:},1) + fit1(vars{:},2)*exp(-x/fit1(vars{:},3));
            exp2 = @(x) fit2(vars{:},1) + fit2(vars{:},2)*exp(-x/fit2(vars{:},3))+ fit2(vars{:},4)*exp(-x/fit2(vars{:},5));
            one_exp_error(force_index,eta_index,alpha_index) = abs(1-exp1(0));
            two_exp_error(force_index,eta_index,alpha_index) = abs(1-exp2(0));
        end
    end
end


for force_index = 1:f
    figure
    surf(X,Y,reshape(one_exp_error(force_index,:,:),[g,a])')
    
    xlabel('eta');
    ylabel('alpha');
    title(['one exponential maximum error, force = ' num2str(forcevec(force_index))])
    set(gca, 'XScale', 'log', 'YScale', 'log');
end


for force_index = 1:f
    figure
    surf(X,Y,reshape(two_exp_error(force_index,:,:),[g,a])')
    
    xlabel('eta');
    ylabel('alpha');
    title(['two exponential maximum error, force = ' num2str(forcevec(force_index))])
    set(gca, 'XScale', 'log', 'YScale', 'log');
end
%%
for force_index = 1:f
    figure
    surf(X,Y,reshape(2*two_exp_error(force_index,:,:),[g,a])'-reshape(one_exp_error(force_index,:,:),[g,a])')
    
    xlabel('eta');
    ylabel('alpha');
    title(['Positive regions are one exponential, force = ' num2str(forcevec(force_index))])
    set(gca, 'XScale', 'log', 'YScale', 'log');
end


%%
for force_index = 1:f
    figure
    surf(X,Y,reshape(error1(force_index,:,:),[g,a])')
    
    xlabel('eta');
    ylabel('alpha');
    title(['one exponential L2 error, force = ' num2str(forcevec(force_index))])
    set(gca, 'XScale', 'log', 'YScale', 'log');
end


for force_index = 1:f
    figure
    surf(X,Y,reshape(error2(force_index,:,:),[g,a])')
    
    xlabel('eta');
    ylabel('alpha');
    title(['two exponential L2 error, force = ' num2str(forcevec(force_index))])
    set(gca, 'XScale', 'log', 'YScale', 'log');
end
%%
error1(:,:,end) = nan*ones(f,g);
error2(:,:,end) = nan*ones(f,g);
 capped_error2 = error2;
% capped_error2(capped_error2<4*1e-5) = 4*1e-5
for force_index = 1:f
%     figure
%     surf(X,Y,reshape(capped_error2(force_index,:,:),[g,a])')
%          
%     xlabel('eta');
%     ylabel('alpha');
%     title(['capped error2, ramptime = ' num2str(forcevec(force_index))])
%      set(gca, 'XScale', 'log', 'YScale', 'log');
    figure
    hold on
    surf(X,Y,reshape(error1(force_index,:,:),[g,a])'./reshape(capped_error2(force_index,:,:),[g,a])')
    alpha(0.7)

    contour(X,Y,reshape(error1(force_index,:,:),[g,a])'./reshape(capped_error2(force_index,:,:),[g,a])',linspace(0,50,11),'ShowText','on')
%     
    xlabel('eta');
    ylabel('alpha');
    title(['L2 error 1 /L2 error 2, ramptime = ' num2str(forcevec(force_index))])
     set(gca, 'XScale', 'log', 'YScale', 'log');
end

%%
error_threshold = 5e-4      ;
colourvec = ['y','m','c','r','b','g','k','y','m','c','r'];%need to stop reuse of colours, temp measure
%colourvec = [linspace(0,1,f);linspace(1,0,f);zeros(1,f)];
figure
hold on
[X,Y] = meshgrid(etavec,alphavec);
contours = logspace(-10,10,21);
plot(X,Y,'ko');
h = [];
for time_index = 1:f
    [~,h(time_index)] = contour(X,Y,reshape(error1(time_index,:,:),[g,a])',[error_threshold,error_threshold],colourvec(time_index),'ShowText','off');
    set(gca, 'XScale', 'log', 'YScale', 'log');
end
xlabel('eta');
ylabel('alpha');
title(['Error contour threshold, ' num2str(error_threshold)])
legendflex(h,arrayfun(@num2str,forcevec,'UniformOutput',false));
%SaveAsPngEpsAndFig(-1,[pwd '/pictures/strainfitting/fitcontour-' num2str(T)]  , 7, 7/5, 9)

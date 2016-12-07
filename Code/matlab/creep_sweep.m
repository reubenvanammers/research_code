%Varies parameters of stress_2d_ode, and plots time strain graphs in a
%subplot.


forcevec = logspace(-1,0,5);
%Tvec = [0 20 100];
T = 50;
gammavec = logspace(-2,0,10);
alphavec = logspace(-2,0,10);
tend = 2000;
f = length(forcevec);g = length(gammavec);a = length(alphavec);
L = f*g*a;
straincell = cell(1,L);
timecell = cell(1,L);
flagcell = cell(1,L);
%restoringcell = cell(1,L);
parfor_progress(L);
parfor i = 1:L%index loops over alpha, then gamma, then T
    counter = i-1;
    alpha_index = mod(counter,a);
    counter = (counter-alpha_index)/a;
    gamma_index = mod(counter,g);
    counter = (counter-gamma_index)/g;
    force_index = counter;
    alpha = alphavec(alpha_index+1);
    force = forcevec(force_index+1);
    gamma = gammavec(gamma_index+1);%converts linear index to alpha,gamma,T
    [Time, Y,~,flag] =stress_2d_ode_maxstrain(alpha,gamma,T,tend,[10,10],force);
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
error_fit = nan*ones(f,g,a);

save([pwd '\workspaces\creeperror']);
%%
load([pwd '\workspaces\creeperror']);
for i = 1:L
    counter = i-1;
    alpha_index = mod(counter,a);
    counter = (counter-alpha_index)/a;
    gamma_index = mod(counter,g);
    counter = (counter-gamma_index)/g;
    force_index = counter;
    vars = {force_index+1,gamma_index+1,alpha_index+1};
    straincell2{vars{:}} = straincell{i};
    timecell2{vars{:}} = timecell{i};
    if ~flagcell{i}
        try
            [fit1(vars{:},:),error1(vars{:}),fit2(vars{:},:),error2(vars{:}),error_fit(vars{:})] = CalculateExponentialFits(timecell2{vars{:}}',straincell2{vars{:}}');
        end
    end
end


%% plots strain-time graphs and exponential fits

for force_index = 1:f
    figure
    hold on
    for alpha_index = 1:a 
        for  gamma_index = 1:g
            subplot(g,a,alpha_index+a*(gamma_index-1))
            vars = {force_index,gamma_index,alpha_index};
            exp1 = @(x) fit1(vars{:},1) + fit1(vars{:},2)*exp(-x/fit1(vars{:},3));
            exp2 = @(x) fit2(vars{:},1) + fit2(vars{:},2)*exp(-x/fit2(vars{:},3))+ fit2(vars{:},4)*exp(-x/fit2(vars{:},5));
            plot(timecell2{vars{:}},straincell2{vars{:}},'r',timecell2{vars{:}},exp1(timecell2{vars{:}}),'k--',timecell2{vars{:}},exp2(timecell2{vars{:}}),'b--')
            title(['alpha = ', num2str(alphavec(alpha_index)), ' gamma = ', num2str(gammavec(gamma_index)), ' force = ' num2str(gammavec(force_index))]); 
        end
    end
end

%%
error_threshold = 1e-5;
colourvec = ['y','m','c','r','b','g','k'];
figure
hold on
[X,Y] = meshgrid(gammavec,alphavec);
contours = logspace(-10,10,21);
plot(X,Y,'ko');

for force_index = 1:f
    [~,h(force_index)] = contour(X,Y,reshape(error_fit(force_index,:,:),[g,a])',[error_threshold,error_threshold],colourvec(force_index),'ShowText','on');
    set(gca, 'XScale', 'log', 'YScale', 'log');
end
xlabel('gamma');
ylabel('alpha');
title('regions above/right of contours are good fits')
legend(h,arrayfun(@num2str,forcevec,'UniformOutput',false));
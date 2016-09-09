Tvec = [0 20 100];
etavec = [1 0.5 0.1];
alphavec = 0.2:0.2:1;
tend = 300;
t = length(Tvec);g = length(etavec);a = length(alphavec);
L = t*a*g;
straincell = cell(1,L);
timecell = cell(1,L);
restoringcell = cell(1,L);
[t,s,r] = strain_restoring(@vertex_ode_reference,1,1,0,1,0,0,300);
parfor_progress(L);
for i = 1:L%index loops over alpha, then eta, then T
    datetime('now')
    counter = i-1
    a1 = mod(counter,a);
    counter = (counter-a1)/a;
    g1 = mod(counter,g);
    counter = (counter-g1)/g;
    t1 = counter;
    alpha = alphavec(a1+1)
    T = Tvec(t1+1)
    eta = etavec(g1+1)%converts linear index to alpha,eta,T
    if alpha == 1;
        straincell{i} = s;
        timecell{i} = t;
        restoringcell{i} = r;
    else
        [Time, strain, restoring] =strain_restoring(@vertex_ode_reference,1,2,0,alpha,eta,T,tend);
        straincell{i} = strain;
        timecell{i} = Time;
        restoringcell{i} = restoring;
    end
    parfor_progress;
end
parfor_progress(0);
strain_matrix = cell2mat(straincell);
time_matrix = cell2mat(timecell);
restoring_matrix = cell2mat(restoringcell);
strainmax = max(max(strain_matrix));
figure %plots strain-time graph and restoringforce-time graphs

for T = 1:t
    for eta = 1:g
        subplot(g,t,T+(t)*(eta-1))
        for alpha = 1:a
            hold on
            i = alpha+(eta-1)*a+(T-1)*a*g;
            h(alpha) = plot(time_matrix(:,i),strain_matrix(:,i),'DisplayName',num2str(alphavec(alpha)));
        end
        axis([0 tend 1 (strainmax+0.1)])
        title([' eta = ',num2str(etavec(eta)),' T = ',num2str(Tvec(T))])
        legend(h)
    end
end
figure            
for T = 1:t
    for eta = 1:g
        subplot(g,t,T+(t)*(eta-1))
        for alpha = 1:a
            hold on
            i = alpha+(eta-1)*a+(T-1)*a*g;
            h(alpha) = plot(time_matrix(:,i),restoring_matrix(:,i),'DisplayName',num2str(alphavec(alpha)));
        end
        title([' eta = ',num2str(etavec(eta)),' T = ',num2str(Tvec(T))])
        legend(h)
    end
end        
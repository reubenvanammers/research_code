gamma = 0;
lambda =3;
eta =1; 
T = 0;
beta=3;
tend = 1000;
i = 0;

alphavec = 0.2:0.2:1;
stressvec = 0.05:0.05:1;
a = length(alphavec); s = length(stressvec);
L = a*s;
straincell = cell(1,L);
stresscell = cell(1,L);
flagcell = cell(1,L);
historycell = cell(1,L);
t_historycell = cell(1,L);
%global external_force monoflag
parfor_progress(L);
parfor i = 1:L;
    counter = i-1;
    s1 = mod(counter,s);
    counter = (counter-s1)/s;
    a1 = mod(counter,a);
    alpha = alphavec(a1+1);
    external_force = stressvec(s1+1);
    [Time, Y,cell_history,cell_t_history,monoflag] = vertex_restructuring_static(lambda,beta,gamma,alpha,eta,T,tend,external_force);
    N = size(Y,2)/4;
    xvalues = Y(:,1:N);
    strain = (max(xvalues,[],2)-min(xvalues(1,:)))/(max(xvalues(1,:))-min(xvalues(1,:)));
    strainmax = strain(end);
    straincell{i} = strainmax;
    stresscell{i} = external_force;
    flagcell{i} = monoflag;
    historycell{i} = cell_history;
    t_historycell{i} =cell_t_history;
    parfor_progress;
end
parfor_progress(0);

alphastrain = cell(1,a);
alphastress = cell(1,a);
alphabreakstrain = cell(1,a);
alphabreakstress= cell(1,a);


for i = 1:L
    counter = i-1;
    s1 = mod(counter,s);
    counter = (counter-s1)/s;
    alpha = mod(counter,a)+1;
    if isempty(alphabreakstrain{alpha})
        alphastrain{alpha} = [alphastrain{alpha} straincell{i}];
        alphastress{alpha} = [alphastress{alpha} stresscell{i}];
    end
    if strcmp(flagcell{i},'break') && isempty(alphabreakstrain{alpha})
        alphabreakstrain{alpha} = [alphabreakstrain{alpha} straincell{i}];
        alphabreakstress{alpha} = [alphabreakstress{alpha} stresscell{i}];
    end
end
hold on
for alpha = 1:a
    h(alpha) = plot(alphastrain{alpha},alphastress{alpha},'DisplayName',num2str(alphavec(alpha)));
    plot(alphabreakstrain{alpha},alphabreakstress{alpha},'kx')
end
legend(h);
title('Stress-Strain t1swap')
xlabel('strain');
ylabel('stress');
% for alpha = 0.2:0.2:1
%     i = i+1;
%     stress_strain = [];
%     for external_force = 0:0.05:0.2
%         [Time, strain, restoring] = strain_restoring(@vertex_restructuring,lambda,beta,gamma,alpha,eta,T,tend);
%         strainmax = strain(end);
%         stress_strain = [stress_strain; external_force strainmax];
%     end
%     
%    stress_strain_cell{i} =  stress_strain;
% end
% 
% figure
% hold on
% for j = 1:length(stress_strain_cell);
%     plot(stress_strain_cell{j}(:,2),stress_strain_cell{j}(:,1));
% end
function [Time,strain,Y] = strain_calc_expfit(fnhandle,varargin)
%calculatse strain for fnhandle with arguments varargin.
[Time,Y]=fnhandle(varargin{:});
N = size(Y,2)/4;
xvalues = Y(:,1:N);
strain = (max(xvalues,[],2)-min(xvalues(1,:)))/(max(xvalues(1,:))-min(xvalues(1,:)));
[fit1,error1,fit2,error2,error_fit_L2,error_fit_inf,~,stretched_params] = CalculateExponentialFits(Time',strain',1,3);

if fit2(3) > fit2(5)
    temptimecoef = fit2(3);
    tempcoef = fit2(2);
    fit2(3) = fit2(5);
    fit2(5) = temptimecoef;
    fit2(2) = fit2(4);
    fit2(4) = tempcoef;
end


exp1 = @(x) fit1(1) + fit1(2)*exp(-x/fit1(3));
exp2 = @(x) fit2(1) + fit2(2)*exp(-x/fit2(3))+ fit2(4)*exp(-x/fit2(5));
exp2_1 = @(x) fit2(1)+fit2(4)+fit2(2)*exp(-x/fit2(3));
exp2_2 = @(x) fit2(1) + fit2(4)*exp(-x/fit2(5));
stretched = @(x) stretched_params(1) + stretched_params(2)*exp(-x.^(stretched_params(4))./stretched_params(3));

figure
plot(Time,strain,'r')

figure
plot(Time,strain,'r',Time,exp1(Time),'k--')
xlabel('Time')
ylabel('Strain')
figure
plot(Time,strain,'r',Time,exp1(Time),'k',Time,exp2(Time),'b--')


figure
plot(Time,strain,'r',Time,exp1(Time),'k',Time,exp2(Time),'b--',Time,exp2_1(Time),'g--',Time,exp2_2(Time),'m--')


figure
plot(Time,strain,'r',Time,stretched(Time),'g--')
stretched_params


figure
global restoring_rec t_rec
l = length(t_rec)-length(restoring_rec);
t_rec = t_rec(l+1:end);
restoring_temp = restoring_rec;
[t_rec,ia,~] = unique(t_rec);
restoring_temp = restoring_temp(ia);%deletes duplicate time entries for interpolation
restoring = interp1(t_rec,restoring_temp,Time);%interpolates restoring force to be same size as time vector
plot(Time(2:end),restoring(2:end))


figure
plot(Time,strain-exp1(Time),'k-',Time,strain-exp2(Time),'b',Time,zeros(size(Time)),'r--')%,Time,strain-exp3(Time),'g-.','MarkerSize',4)            %plot(timecell3{vars{1:end-1}},straincell3{vars{1:end-1}},'r',timecell3{vars{1:end-1}},exp1(timecell3{vars{1:end-1}}),'k--',timecell3{vars{1:end-1}},exp2(timecell3{vars{1:end-1}}),'b--',timecell3{vars{1:end-1}},exp2_1(timecell3{vars{1:end-1}}),'g-.',timecell3{vars{1:end-1}},exp2_2(timecell3{vars{1:end-1}}),'y-.')
xlabel('Time')
ylabel('Strain Approx Difference')



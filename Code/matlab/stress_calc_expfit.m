%calculates stress for fnhandle with arguments varargin, and then fits and
%plots with various potential fits, such as one, two, stretched
%exponentials etc.

function [Time,stress,Y] = stress_calc_expfit(fnhandle,varargin)
[Time,Y,~,stress,t_rec,stress_index]=fnhandle(varargin{:});

stress



tstop = t_rec(stress_index)
Time = Time-tstop
Time = Time(Time >=0)

stress = stress(stress_index:end);
t_rec = t_rec(stress_index:end)-tstop;

stress_temp = stress;
[t_rec,ia,~] = unique(t_rec);


stress_temp = stress_temp(ia)%deletes duplicate time entries for interpolation


stress2 = interp1(t_rec,stress_temp,Time);%interpolates stress force to be same size as time vector

size(Time)
size(stress2)


[fit1,error1,fit2,error2,error_fit_L2,error_fit_inf,~,stretched_params] = CalculateExponentialFits(Time',stress2',-1,1);

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
plot(Time,stress2,'r')

figure
plot(Time,stress2,'r',Time,exp1(Time),'k--')
xlabel('Time')
ylabel('stress')
figure
plot(Time,stress2,'r',Time,exp1(Time),'k',Time,exp2(Time),'b--')


figure
plot(Time,stress2,'r',Time,exp1(Time),'k',Time,exp2(Time),'b--',Time,exp2_1(Time),'g--',Time,exp2_2(Time),'m--')


figure
plot(Time,stress2,'r',Time,stretched(Time),'g--')
stretched_params


figure
plot(Time,stress2-exp1(Time),'k-',Time,stress2-exp2(Time),'b',Time,zeros(size(Time)),'r--')%,Time,stress-exp3(Time),'g-.','MarkerSize',4)            %plot(timecell3{vars{1:end-1}},stresscell3{vars{1:end-1}},'r',timecell3{vars{1:end-1}},exp1(timecell3{vars{1:end-1}}),'k--',timecell3{vars{1:end-1}},exp2(timecell3{vars{1:end-1}}),'b--',timecell3{vars{1:end-1}},exp2_1(timecell3{vars{1:end-1}}),'g-.',timecell3{vars{1:end-1}},exp2_2(timecell3{vars{1:end-1}}),'y-.')
xlabel('Time')
ylabel('stress Approx Difference')


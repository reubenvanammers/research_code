function [Time,strain,Y] = strain_calc_expfit(fnhandle,varargin)
%calculatse strain for fnhandle with arguments varargin.
[Time,Y]=fnhandle(varargin{:});
N = size(Y,2)/4;
xvalues = Y(:,1:N);
strain = (max(xvalues,[],2)-min(xvalues(1,:)))/(max(xvalues(1,:))-min(xvalues(1,:)));
[fit1,error1,fit2,error2,error_fit_L2,error_fit_inf,~,stretched_params] = CalculateExponentialFits(Time',strain',1,3);

exp1 = @(x) fit1(1) + fit1(2)*exp(-x/fit1(3));
exp2 = @(x) fit2(1) + fit2(2)*exp(-x/fit2(3))+ fit2(4)*exp(-x/fit2(5));
exp2_1 = @(x) fit2(1)+fit2(4)+fit2(2)*exp(-x/fit2(3));
exp2_2 = @(x) fit2(1) + fit2(4)*exp(-x/fit2(5));
stretched = @(x) stretched_params(1) + stretched_params(2)*exp(-x.^(stretched_params(4))./stretched_params(3));

figure
plot(Time,strain,'r')

figure
plot(Time,strain,'r',Time,exp1(Time),'k',Time,exp2(Time),'b-')


figure
plot(Time,strain,'r',Time,exp1(Time),'k',Time,exp2(Time),'b-',Time,exp2_1(Time),'g--',Time,exp2_2(Time),'m--')


figure
plot(Time,strain,'r',Time,stretched(Time),'g--')
stretched_params

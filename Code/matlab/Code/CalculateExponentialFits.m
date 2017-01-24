function [ one_exp_params, one_exp_ave_error ,two_exp_params, two_exp_ave_error, error_fit] = CalculateExponentialFits(time,data)
%CalculateExponentialFits Calclate the fit and error to exponential models

options_1 = fitoptions('Method', 'NonLinearLeastSquares','TolFun',1e-9)%,'Algorithm','Levenberg-Marquardt');

% This is an exact fit and used as the initial guess for the
% iterative solver
temp_fit = exp2fit(time,data,1);

initial = [temp_fit(1),temp_fit(2),temp_fit(3)];
options_1 = fitoptions('Method', 'NonLinearLeastSquares','Algorithm','Levenberg-Marquardt','Start',initial,'TolFun',1e-9);
oneexpfit = fittype('a + b*exp(-x/c)','dependent',{'y'},'independent',{'x'},'coefficients',{'a','b','c'})
[f,error] = fit( time', data', oneexpfit,options_1);
one_exp_params = [f.a,f.b,f.c];
one_exp_ave_error = error.sse/length(time);

initial = [temp_fit(1),temp_fit(2),temp_fit(3),1,1]
options_2 = fitoptions('Method', 'NonLinearLeastSquares','Algorithm','Levenberg-Marquardt','Start',initial,'TolFun',1e-9);
twoexpfit = fittype('a + b*exp(-x/c) + d*exp(-x/e)','dependent',{'y'},'independent',{'x'},'coefficients',{'a','b','c','d','e'})
[f,error] = fit( time', data', twoexpfit,options_2);
two_exp_params = [f.a,f.b,f.c,f.d,f.e];
two_exp_ave_error = error.sse/length(time);


one_exp_curve = one_exp_params(1)+one_exp_params(2).*exp(-time./one_exp_params(3));
two_exp_curve = two_exp_params(1)+two_exp_params(2).*exp(-time./two_exp_params(3))+two_exp_params(4).*exp(-time./two_exp_params(5));

error_fit = sum((one_exp_curve-two_exp_curve).^2)/length(time);

end


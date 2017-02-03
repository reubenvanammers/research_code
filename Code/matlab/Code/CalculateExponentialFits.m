function [ one_exp_params, one_exp_ave_error ,two_exp_params, two_exp_ave_error, error_fit,error_integrate] = CalculateExponentialFits(time,data)
%CalculateExponentialFits Calclate the fit and error to exponential models

options_1 = fitoptions('Method', 'NonLinearLeastSquares','TolFun',1e-9)%,'Algorithm','Trust-Region');

% This is an exact fit and used as the initial guess for the
% iterative solver
temp_fit = exp2fit(time,data,1);

initial = [temp_fit(1),temp_fit(2),temp_fit(3)];
options_1 = fitoptions('Method', 'NonLinearLeastSquares','Algorithm','Trust-Region','Start',initial,'TolFun',1e-9);

oneexpfit = fittype('a + b*exp(-x/c)','dependent',{'y'},'independent',{'x'},'coefficients',{'a','b','c'})
[f,error] = fit( time', data', oneexpfit,options_1);
one_exp_params = [f.a,f.b,f.c];
one_exp_ave_error = error.sse/length(time);

initial = [temp_fit(1),temp_fit(2),temp_fit(3),1,1]
%options_2 = fitoptions('Method', 'NonLinearLeastSquares','Algorithm','Trust-Region','Start',initial,'TolFun',1e-9);
options_2 = fitoptions('Method', 'NonLinearLeastSquares','Algorithm','Trust-Region','Start',initial,'TolFun',1e-9,'Lower',[-Inf -Inf 0 -Inf 0]);

twoexpfit = fittype('a + b*exp(-x/c) + d*exp(-x/e)','dependent',{'y'},'independent',{'x'},'coefficients',{'a','b','c','d','e'})
[f,error] = fit( time', data', twoexpfit,options_2);
two_exp_params = [f.a,f.b,f.c,f.d,f.e];
two_exp_ave_error = error.sse/length(time);


one_exp_curve = one_exp_params(1)+one_exp_params(2).*exp(-time./one_exp_params(3));
two_exp_curve = two_exp_params(1)+two_exp_params(2).*exp(-time./two_exp_params(3))+two_exp_params(4).*exp(-time./two_exp_params(5));

error_fit = sum((one_exp_curve-two_exp_curve).^2)/length(time);


C1 = one_exp_params(1);A1 = one_exp_params(2); B1 = one_exp_params(3); C2 = two_exp_params(1); A2 = two_exp_params(2); B2 = two_exp_params(3); D2 =two_exp_params(4); G2 = two_exp_params(5);

% integral of (f-g)^2 over the region 0 and 1, where f is one_exp_curve and
% g is two_exp_curve
error_integrate = (1/2).*(A1.^2.*B1+A2.^2.*B2+4.*A1.*B1.*(C1+(-1).*C2)+2.*(C1+(-1).* ...
  C2).^2+(-4).*A2.*B2.*(A1.*B1.*(B1+B2).^(-1)+C1+(-1).*C2)+(-1).* ...
  A1.^2.*B1.*exp(1).^((-2).*B1.^(-1))+(-4).*A1.*B1.*(C1+(-1).*C2).* ...
  exp(1).^((-1).*B1.^(-1))+(-1).*A2.^2.*B2.*exp(1).^((-2).*B2.^(-1)) ...
  +4.*A2.*B2.*exp(1).^((-1).*B2.^(-1)).*(C1+(-1).*C2+A1.*B1.*(B1+B2) ...
  .^(-1).*exp(1).^((-1).*B1.^(-1)))+D2.^2.*G2+(-1).*D2.^2.*exp(1).^( ...
  (-2).*G2.^(-1)).*G2+(-4).*D2.*G2.*(C1+(-1).*C2+A1.*B1.*(B1+G2).^( ...
  -1)+(-1).*A2.*B2.*(B2+G2).^(-1))+4.*D2.*exp(1).^((-1).*G2.^(-1)).* ...
  G2.*(C1+(-1).*C2+A1.*B1.*exp(1).^((-1).*B1.^(-1)).*(B1+G2).^(-1)+( ...
  -1).*A2.*B2.*exp(1).^((-1).*B2.^(-1)).*(B2+G2).^(-1)));

%error_integrate = max(abs(one_exp_curve-two_exp_curve));

end


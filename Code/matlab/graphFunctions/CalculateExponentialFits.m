function [ one_exp_params, one_exp_ave_error ,two_exp_params, two_exp_ave_error, error_fit,infinity_error,three_exp_params,stretched_params] = CalculateExponentialFits(time,data,exponential_sign,guesstype)
%CalculateExponentialFits Calclate the fit and error to exponential models
%exponential sign is posiitve if fits are such that exponentials should be
%increasing, 0 either way, and -1 if exponential should be decreasing
%guesstype determines the initial condition used for the guess
if nargin <4
    guesstype = 1;
end

if nargin <3
    exponential_sign = 0;
end

if exponential_sign == -1
    lower_bounds = [-Inf 0 0 0 0];
    upper_bounds = [Inf Inf Inf Inf Inf];
elseif exponential_sign ==0
    lower_bounds = [-Inf -Inf 0 -Inf 0];
    upper_bounds = [Inf Inf Inf Inf Inf];
elseif exponential_sign == 1
    lower_bounds = [-Inf -Inf 0 -Inf 0];
    upper_bounds = [Inf 0 Inf 0 Inf];
end

    

options_1 = fitoptions('Method', 'NonLinearLeastSquares','TolFun',1e-9)%,'Algorithm','Trust-Region');

% This is an exact fit and used as the initial guess for the
% iterative solver
temp_fit = exp2fit(time,data,1);

if guesstype ==1 || guesstype ==2
    initial_1 = [temp_fit(1),temp_fit(2),temp_fit(3)];
else
    initial_1 = [0, -1, 1];
end

%initial_1 = [temp_fit(1),temp_fit(2),temp_fit(3)];
%options_1 = fitoptions('Method', 'NonLinearLeastSquares','Algorithm','Levenberg-Marquardt','Start',initial,'TolFun',1e-9);
options_1 = fitoptions('Method', 'NonLinearLeastSquares','Algorithm','Trust-Region','Start',initial_1,'TolFun',1e-9,'Lower',[-Inf -Inf 0],'Upper',[Inf Inf Inf]);

oneexpfit = fittype('a + b*exp(-x/c)','dependent',{'y'},'independent',{'x'},'coefficients',{'a','b','c'})
[f,error] = fit( time', data', oneexpfit,options_1);
one_exp_params = [f.a,f.b,f.c];
one_exp_ave_error = error.sse/length(time);

if guesstype ==1
    initial_2 = [temp_fit(1),temp_fit(2),temp_fit(3),1,1];
    % Initial guess has one exponential fti with peturbation away from it;
elseif guesstype ==2
    initial_2 = [temp_fit(1),temp_fit(2)/2,temp_fit(3),temp_fit(2)/2,temp_fit(3)];
    %two exp initial guess is one exponential guess, each exponential
    %equally weighted
elseif guesstype ==3
    %initial_2 = [0 1 1 1 1];    % Uses guess not related to 1 exponential guess
    if exponential_sign == 1
        initial_2 = [rand -rand rand -rand rand];
    else
       initial_2 = rand(1,5);
    end
end
    

%initial = [temp_fit(1),temp_fit(2)/2,temp_fit(3),temp_fit(2)/2,temp_fit(3)]
%initial_2 = [temp_fit(1),temp_fit(2),temp_fit(3),1,1]

%options_2 = fitoptions('Method', 'NonLinearLeastSquares','Algorithm','Levenberg-Marquardt','Start',initial,'TolFun',1e-9);
options_2 = fitoptions('Method', 'NonLinearLeastSquares','Algorithm','Trust-Region','Start',initial_2,'TolFun',1e-9,'Lower',lower_bounds,'Upper',upper_bounds)

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
% error_fit = (1/2).*(A1.^2.*B1+A2.^2.*B2+4.*A1.*B1.*(C1+(-1).*C2)+2.*(C1+(-1).* ...
%   C2).^2+(-4).*A2.*B2.*(A1.*B1.*(B1+B2).^(-1)+C1+(-1).*C2)+(-1).* ...
%   A1.^2.*B1.*exp(1).^((-2).*B1.^(-1))+(-4).*A1.*B1.*(C1+(-1).*C2).* ...
%   exp(1).^((-1).*B1.^(-1))+(-1).*A2.^2.*B2.*exp(1).^((-2).*B2.^(-1)) ...
%   +4.*A2.*B2.*exp(1).^((-1).*B2.^(-1)).*(C1+(-1).*C2+A1.*B1.*(B1+B2) ...
%   .^(-1).*exp(1).^((-1).*B1.^(-1)))+D2.^2.*G2+(-1).*D2.^2.*exp(1).^( ...
%   (-2).*G2.^(-1)).*G2+(-4).*D2.*G2.*(C1+(-1).*C2+A1.*B1.*(B1+G2).^( ...
%   -1)+(-1).*A2.*B2.*(B2+G2).^(-1))+4.*D2.*exp(1).^((-1).*G2.^(-1)).* ...
%   G2.*(C1+(-1).*C2+A1.*B1.*exp(1).^((-1).*B1.^(-1)).*(B1+G2).^(-1)+( ...
%   -1).*A2.*B2.*exp(1).^((-1).*B2.^(-1)).*(B2+G2).^(-1)));

infinity_error = max(abs(one_exp_curve-two_exp_curve));

if exponential_sign == -1
    lower_bounds = [-Inf 0 0 0 0 0 0 ];
    upper_bounds = [Inf Inf Inf Inf Inf Inf Inf];
elseif exponential_sign ==0
    lower_bounds = [-Inf -Inf 0 -Inf 0 -Inf 0];
    upper_bounds = [Inf Inf Inf Inf Inf Inf Inf];
elseif exponential_sign == 1
    lower_bounds = [-Inf -Inf 0 -Inf 0 -Inf 0];
    upper_bounds = [Inf 0 Inf 0 Inf 0 Inf];
end
initial_3 = [rand,-rand,rand,-rand,rand,-rand,rand];
options_3 = fitoptions('Method', 'NonLinearLeastSquares','Algorithm','Trust-Region','Start',initial_3,'TolFun',1e-9,'Lower',lower_bounds,'Upper',upper_bounds);

threeexpfit = fittype('a + b*exp(-x/c) + d*exp(-x/e)+f*exp(-x/g)','dependent',{'y'},'independent',{'x'},'coefficients',{'a','b','c','d','e','f','g'});
[fit_data,error] = fit( time', data', threeexpfit,options_3);
three_exp_params = [fit_data.a,fit_data.b,fit_data.c,fit_data.d,fit_data.e,fit_data.f,fit_data.g];
three_exp_ave_error = error.sse/length(time);
three_exp_curve = three_exp_params(1)+three_exp_params(2).*exp(-time./three_exp_params(3))+three_exp_params(4).*exp(-time./three_exp_params(5))+three_exp_params(6).*exp(-time./three_exp_params(7));

%stretched exponential


stretchedoptions = fitoptions('Method', 'NonLinearLeastSquares','Algorithm','Trust-Region','Start',[0 -1 1 1],'TolFun',1e-9,'Lower',[-Inf -Inf 0 0],'Upper',[Inf Inf Inf 5]);

stretchedfit = fittype('a + b*exp(-x^p/c)','dependent',{'y'},'independent',{'x'},'coefficients',{'a','b','c','p'});
[fit_data,error] = fit( time', data',stretchedfit,stretchedoptions);
stretched_params = [fit_data.a,fit_data.b,fit_data.c,fit_data.p];


end

    
function linspring_stress_strain(alpha,eta,T,tend,strainfunc,varargin)
global t_rec stress_rec

[Time,Y,Tri] = strain_2d_ode(alpha,eta,T,tend,strainfunc,varargin{:});
l = length(t_rec)-length(stress_rec);
t_rec2 = t_rec(l+1:end);
plot(strainfunc(t_rec2),stress_rec)
xlabel('Strain')
ylabel('Stress')


% figure
% N = size(Y,2)/4;
% xvalues = Y(:,1:N);
% strain = (max(xvalues,[],2)-min(xvalues(1,:)))/(max(xvalues(1,:))-min(xvalues(1,:)));
% plot(Time,strain,Time,strainfunc(Time));
% norm(strain'-strainfunc(Time),Inf)
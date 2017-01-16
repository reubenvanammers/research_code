strain_min = 1.2;
strain_max = 1.6;
tend = 200;
cycle_length = 20;
alpha = 0.5;
eta =0.1;
T = 20;
strainfunc = ramp2(strain_min,strain_max,inf,cycle_length);
[Time,Y,Tri2,stress_rec2]=strain_2d_ode(alpha,eta,T,tend,strainfunc,inf,[10 10]);   
figure
plot(Time,stress_rec2,'k',Time,strainfunc(Time),'r')
legend(['Stress';'Strain'])
xlabel('Time')
figure
plot(strainfunc(Time),stress_rec2)
xlabel('Strain')
ylabel('Stress')
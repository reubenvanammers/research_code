%% Cellular Remodelling Code
% This Matlab code repository is to run dynamic discrete cell simulations
% with remodelling via the use of a reference state. This is implemented in
% various ways for both cell centre and vertex based modelsl. 


%% Cell Centre Models
%Cellular Remodelling has been implemented via cell centre linear spring
%models. Code has been written to for both applying constant stress to the
%righthand side of a monolayer while keeping the left hand fixed, and also
%forcing the monolayer to obey a given stress function. 


% To run a constant stress simulation:
%[Time,Y,Tri,flag]=stress_2d_ode(alpha,eta,T_mem,t_end,gridsize,ext_force,maxstrain)
%Input Variables:
%alpha and eta are parameters of the model, and details are given in the
%paper. T_mem gives the averaging time, with T_mem = 0 for no memory. t_end
%gives the final time of the simulation. 
%Optional Arguments: gridize is a 1*2 vector specifying the dimensions of 
%the monolayer, ext_force is the force applied to the RHS of the monolayer,
%defaulting to 0.2, and maxstrain is a argument to stop the simulation if
%the strain goes above the given value. 

%Output arguments: Time is the output times of the solver, Y is the
%location of the set of points for each time in T for both the real and
%reference states, Tri specifies the mesh of the triangles, and flag is an
%argument giving true if the solver stops prematurely via maxstrain.

%The simulation for both the real and reference states can be visualized
%via the function:
%tri_vis(Time,Y,Tri). eg:
[Time,Y,Tri] = stress_2d_ode(0.5,0.5,0,100);
tri_vis(Time,Y,Tri);

%The strain can also be visualised by passing the function into
%strain_calc with the arguments to be passed in:
%[Time,strain,Y] = strain_calc(fnhandle,varargin)
%eg 
[Time,strain] = strain_calc(@stress_2d_ode,0.5,0.5,0,100);
plot(Time,strain);
%this strain_calc function can take any appropriate function to calculate
%the strains, including for the vertex based solvers. 
%%

%The strain can also be forced on the righthand side of the monolayer with
%a similar function:
%[Time,Y,Tri2,stress_rec] = strain_2d_ode(alpha,eta,T_mem,t_end,strainfunc,t_strain_end,gridsize)
%the arguments are similar to the previous function, but with strainfunc
%giving the function of the strain to force the monolayer into, and
%optional argument t_strain_end to allow the monolayer tor relax after this
%time. Stress_rec records stress experienced on the righthand side of the
%monolayer following this loading. 
%For the stress relaxation experiment, the required force function has been
%implemented via the ramp function ramp(strainmax,reps,reptime), with
%strainmax the maximum strain, reps the number of repetitions of this
%loading (1 for stress relaxation experiment) and reptime is the time it
%takes to fully move to the maximum strainvalue. Eg

[Time,Y,Tri2]=strain_2d_ode(0.5,0.5,0,100,ramp(1.5,1,20));
tri_vis(Time,Y,Tri);

%The strain increases for 20 time units, at which point the righthand side
%of the real state's monolayer stops moving, allowing the system to relax.

%Cyclic loading experiments can be run by choosing an appropriate strain
%function to input into the model, eg a sinusoidal based function. 

%%
%In order to run sweeps of such simulations in order to view compound plots
%to view the behaviour of such simulations, we can run sweep driver files
%for both the creep and strain relaxation experiments.

%For the creep experiment, this is located in the file
%creep_sweep_origtimes.m file. This will run the the simulations, and save
%the time-strain data for each parameter combination in a high dimensional
%matlab cell. The data is then fitted with a variety of fits (via the
%CalculateExponentialFits function, and then this information is then used
%to create contour and surface plots of various features, such as the 
%coefficient ratio of the data. 

%To choose the bounds of such a contour sweep, simply change the bounds at
%the start of the file

%fbounds = {-1 0 3}; %F_ext value
%Tbounds = {1 2 2}; %T_mem value
%etabounds = {-1 0 21};
%alphabounds = {-1 0 21};

%The first argument of these bounds is the log10 value of the minimum
%value, the second argument is the log10 of maximum value, and the final
%value is the number of points to sample between these (inclusive) bounds
%in an exponential spacing. Such a sweep is computationally expensive, so
%running this sweep on multiple threads is recommended.
%The indices of the cells of the various data is arranged in the order
%(F_ext,eta,alpha,T_mem), with the index value being the location of the
%given value in the above bound values. 

%This can be done in similar manner for the stress relaxation experiment
%for recording the time-stress values, in the file strain_sweep_origtimes.m.
%This has similar lines for choosing bounds, but with fbounds being
%replaced with 
%Tbounds = {0 2 3}, with the same syntax as previously.
%Again, various plots, such as for the coefficient ratio for plots
%regarding the time-stress values for the simulations can then be run
%further down in the file.




%% Vertex based models
% There is also a significant amount of code to run simulations for vertex
% based models, with a similar format to the cell centre based simulations
% for both the stress relaxation and creep experiments, all based on a
%Nagai-Honda model of force. With a reference state that has area and
%circumference coupled, simulations for the creep experiment can be run as
%follows:

%[Time,Y,C,flag] =
%creep_vertex(lambda,beta,gamma,alpha,eta,T_mem,tend,gridsize,ext_force,maxstrain)
%Many fo the arguments are similar to the cell centre model described
%above, but the extra complexity of the vertex model requires additional
%parameters to be used; lambda, beta and gamma.
%(Implementing vertex dynamics models of cell populations in biology within
%a consistent computational framework, Progress in Biophysics and Molecular
%Biology) lambda gives the strenght of the volume target constraint, beta
%the circumference based constraint, and gamma the strength of the
%conribution from cell-cell adhesion energy. For this model with the
%reference state, gamma should generally be set to 0 to have physically
%meaningful interactions with the reference state. 
%The output variable C is information regarding the ownership of vertices
%within cells, needed for visualisation of the simualation output, and Y is
%the locations of the cell vertices at each location. 
%Simulations can then be %visualised via the hex_vis_2 function: 
%hex_vis_2(Time,Y,C).

%For example:
[Time,Y,C] = creep_vertex(100,10,0,0.5,0.5,0,100);
hex_vis_2(Time,Y,C);
%%
%Similarly, the stress relaxation can be run via the use of the
%strain_vertex function. 
%[Time,Y,C,stress_rec] = strain_vertex(lambda,beta,gamma,alpha,eta,T,tend,gridsize,strainfunc,t_ramp_end,t_strain_end2)
eg 
[Time,Y,C] = strain_vertex(1,1,0,0.5,0.5,0,200,ramp(1.2,1,100));
hex_vis_2(Time,Y,C);

%In these vertex based models, very small edges in cell vertices can cause
%ode solvers to have trouble. Therefore, without cell rearrangment such as
%T1 swaps, solvers may stall at high forces/quick loading. 


%% Vertex Model cell directionality
%We can also have cell models where there is an additional directional
%component to the cell model along the major axis, chosen to be in the
%direction of cell movement in the positive x direction. All remodelling
%goes through this component, while area and circumference components have
%a static component. 
%[Time,Y,C] = creep_vertex_axis(lambda,beta,gamma,delta,alpha,eta,T,tend,gridsize,ext_force,maxstrain);
%In this model, gamma now represents the strength of this directional component,
%and delta represents the strength of a rotational componenent whose effects
%that hasn't been fully tested, and so generally set this value to 0. 
%eg
[Time,Y,C2] = creep_vertex_axis(100,10,50,0,0.5,0.5,0,100);
hex_vis_2(Time,Y,C);
%Similarly, strain_vertex_axis represents this directional cell model with
%an external strain function.


%% Cell Restructuring
%There are multiple operations that can be done in order change cell shape
%for the reference state. T1 swaps have been implemented for the reference
%state via polygon vertices being too close in the real state. 

%This is implemented for the creep in experiment in the functions 
%vertex_restructuring and vertex_restructuring_static, the latter 
%outputting a flag if the monolayer reaches equilbrium or breaks into two
%seperate segments. 

%As the cell configuration will change based on these T1 swaps, it
%necessiates outputting the history of the cell configuration. This is done
%via the output cells, cell_history, and the times of each of the cell
%configurations, the vector cell_t_history.
%

%[Time,Y,cell_history,cell_t_history] = vertex_restructuring(lambda,beta,gamma,alpha,eta,T,tend)

%This can then be inputted into  hex_vis_2 to visualise the results

[Time,Y,cell_history,cell_t_history] = vertex_restructuring(1,1,0,1,0,0,5); %find an example which changes cell connectivity
hex_vis_2(Time,Y,cell_history,cell_t_history);
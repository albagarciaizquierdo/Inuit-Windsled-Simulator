clc;
clear all;
close all;

% Define p as global
global p
frame = 0;

% Import parameters
p = parameters;

% Plot settings
set(groot,'defaultAxesFontSize', 11);
set(groot,'defaultLineLineWidth', 1);
set(groot,'defaultAxesColorOrder', [0 0 1;1 0 0; 0 0.6 0;1 0.7 0.1; 0.5 0.2 0.6]);
set(groot,'defaultAxesTickLabelInterpreter', 'latex');
set(groot,'defaultAxesGridAlpha', 0.3);
set(groot,'defaultAxesXMinorTick', 'on');
set(groot,'defaultAxesYMinorTick', 'on');
set(groot,'defaultFigureRenderer', 'painters');
set(groot,'defaultLegendBox', 'off');
set(groot,'defaultLegendInterpreter', 'latex');
set(groot,'defaultLegendLocation', 'best');
set(groot,'defaultLineMarkerSize', 1);
set(groot,'defaultTextInterpreter', 'latex');

%% EQUILIBRIUM STATE CALCULATION
% Initial guesses for equilibrium
Gamma0 = 60*acos(-1)/180;
alpha0 = 10*acos(-1)/180;
x_K0 = - (p.l0+1)*cos(Gamma0);
z_K0 = - (p.l0+1)*sin(Gamma0);
x_red0 = [x_K0 z_K0 alpha0]';

% Equilibirum state calculation
[X_red_eq, ~, ~] = my_fzero("fun_equilibrio_red",x_red0,1e-8,30,1e-6)
x_eq = equilibrium_conditions(X_red_eq);
[gamma_p_eq, gamma_m_eq] = fun_gamma(x_eq,p);

% Equilibrium state drawing
title_eq = "System in equilibrium";
box = 0;
fun_draw_system2(p,x_eq,title_eq,box)

% Equilibrium state results
% [F_S,M_OS,T_ASp,T_ASm,W_S,N,F_r,ASp_AKp,ASm_AKm,up,um,ASp_OS,ASm_OS,OK_AKp,OK_AKm,F_K,M_OK,v_OK,omega_KE,H_OK,W_K,T_AKp,T_AKm,F_a,M_a,alpha,beta,v_a,l_p,l_m,dxdt_eq] = fun_get_results(x_eq,p,frame)

%% STABILITY
% Jacobian, eigenvalues and eigenvectors
dh = 1e-6;
J_eq = fun_jac_num(@RHS,0,x_eq,dh);

imag_J = any(imag(J_eq(:)) ~= 0);
if imag_J
    disp('Jacobian has imaginary values');
else
    disp('Jacobian has NO imaginary values');
end

[vec_J, val_J]=eig(J_eq);

format long
eigenvalues_J = diag(val_J)

J_cond = cond(J_eq)
rho = zeros(1,18);
omega = zeros(1,18);

for i = 1:length(eigenvalues_J)
    rho(i) = real(eigenvalues_J(i));
    omega(i) = imag(eigenvalues_J(i));
end

%% INITIAL CONDITIONS
x0 = initial_conditions;

% draw_system(p,x0)
% title('Initial conditions')

%% INTEGRATOR CHOICE
% steps = 500;
% tspan = linspace(1,100,steps);
% opts = odeset('RelTol',1e-6,'AbsTol',1e-8);
% tic;
% [t,x] = ode45(@RHS,tspan,x_eq,opts);
% time_ode45 = toc;
% 
% [t,x] = ode23(@RHS,tspan,x_eq,opts);
% time_ode23 = toc;
% 
% [t,x] = ode113(@RHS,tspan,x_eq,opts);
% time_ode113 = toc;
% 
% [t,x] = ode89(@RHS,tspan,x_eq,opts);
% time_ode89 = toc;
% 
% tic;
% [t,x] = ode15s(@RHS,tspan,x_eq,opts);
% time_ode15s = toc;
% 
% tic;
% [t,x] = ode23s(@RHS,tspan,x_eq,opts);
% time_ode23s = toc;
% 
% tic;
% [t,x] = ode23t(@RHS,tspan,x_eq,opts);
% time_ode23t = toc;
% 
% tic;
% [t,x] = ode23tb(@RHS,tspan,x_eq,opts);
% time_ode23tb = toc;

%% INTEGRATOR
tic;
steps = 500;
tspan = linspace(1,100,steps);
opts = odeset('RelTol',1e-7,'AbsTol',1e-7);
% [t,x] = ode15s(@RHS,tspan,x_eq,opts);
% [t,x] = ode15s(@RHS,tspan,x_eq+1e-3*rand(18,1),opts);
% [t,x] = ode15s(@RHS,tspan,x_eq,opts);
[t,x] = ode15s(@RHS,tspan,x0,opts);

time_ode = toc;
%% RESULTS
[F_S,M_OS,T_ASp,T_ASm,W_S,N,F_r,v_OS,ASp_AKp,ASm_AKm,up,um,ASp_OS,ASm_OS,OK_AKp,...
    OK_AKm,F_K,M_OK,v_OK,omega_KE,H_OK,W_K,T_AKp,T_AKm,F_a,M_a,alpha,beta,v_a,l_p,l_m,dxdt] = fun_get_results(x,p,frame);

%% PLOTS
fun_draw_results_eq(x,t,p,frame)

%% VIDEO
% fun_video(p,x,"xeq_50_fixed",t)

%% SAVE FIGURES
% % Create folder
% folderName = 'Figures_04_07_pert_ve-4';
% if ~exist(folderName, 'dir')
%     mkdir(folderName);
% end
% 
% % Obtain figures object
% figures = findall(groot, 'Type', 'figure');
% 
% % Save each figure as jpg in the created folder
% for i = 1:9
%     fig = figures(i); 
%     filename = sprintf('%s/figure_%d.jpg', folderName, i);
%     saveas(fig, filename);
% end

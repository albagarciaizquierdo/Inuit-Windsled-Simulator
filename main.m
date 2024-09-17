%% A Simulation Tool for the Dynamic and Stability of the Inuit Windsled
% Alba Garcia Izquierdo
% September 2024

%%
clc;
clear all;
close all;

% Define p as global
global p
frame = 0; % All results to be shown in S_E

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
fprintf('EQUILIBRIUM STATE CALCULATION \n');
% Initial guesses for equilibrium
gamma0 = 60*acos(-1)/180;
alpha0 = 10*acos(-1)/180;
x_K0 = - 1.1*p.l0*cos(gamma0);
z_K0 = - 1.1*p.l0*sin(gamma0);
x_red0 = [x_K0 z_K0 alpha0]';
fprintf('Initial guesses: \n x_eq0 = %.2f m \n z_eq0 = %.2f m \n theta_eq0 = %.4f rad = %.2f deg \n', x_K0, z_K0, alpha0, rad2deg(alpha0));

% Equilibirum state calculation
[X_red_eq, ~, ~] = my_fzero("fun_equilibrio_red",x_red0,1e-8,30,1e-6);
x_eq = fun_equilibrium_conditions(X_red_eq);
[gamma_p_eq, gamma_m_eq] = fun_gamma(x_eq,p);
fprintf('Equilibrium state: \n x_eq = %.2f m \n z_eq = %.2f m \n theta_eq = %.4f rad = %.2f deg \n gamma_eq = %.2f deg \n', X_red_eq(1), X_red_eq(2), X_red_eq(3), rad2deg(X_red_eq(3)), gamma_p_eq);

% Equilibrium state drawing
title_eq = "System in equilibrium";
box = 0;
fun_draw_system2(p,x_eq,title_eq,box)

% Equilibrium state results
% [F_S,M_OS,T_ASp,T_ASm,W_S,N,F_r,ASp_AKp,ASm_AKm,up,um,ASp_OS,ASm_OS,OK_AKp,OK_AKm,F_K,M_OK,v_OK,omega_KE,H_OK,W_K,T_AKp,T_AKm,F_a,M_a,alpha,beta,v_a,l_p,l_m,dxdt_eq] = fun_get_results(x_eq,p,frame)
fprintf('\n');
%% ODE SOLVER SELECTION
fprintf('ODE SOLVER SELECTION \n');

steps = 500;
tspan = linspace(1,100,steps);
opts = odeset('RelTol',1e-6,'AbsTol',1e-8);
x_per = x_eq + 1e-3;

tic;
[t,x] = ode45(@RHS,tspan,x_per,opts);
time_ode45 = toc;

[t,x] = ode23(@RHS,tspan,x_per,opts);
time_ode23 = toc;

[t,x] = ode113(@RHS,tspan,x_per,opts);
time_ode113 = toc;

tic;
[t,x] = ode15s(@RHS,tspan,x_per,opts);
time_ode15s = toc;

tic;
[t,x] = ode23s(@RHS,tspan,x_per,opts);
time_ode23s = toc;

tic;
[t,x] = ode23t(@RHS,tspan,x_per,opts);
time_ode23t = toc;

tic;
[t,x] = ode23tb(@RHS,tspan,x_per,opts);
time_ode23tb = toc;

time_odes = [time_ode45, time_ode23, time_ode113, time_ode15s, time_ode23s, time_ode23t, time_ode23tb];
solvers = {'ode45', 'ode23', 'ode113', 'ode15s', 'ode23s', 'ode23t', 'ode23tb'};
for i = 1:length(solvers)
    fprintf(' %s time: %.2f seconds \n', solvers{i}, time_odes(i));
end
[min_time, idx] = min(time_odes);

fprintf('The fastest ODE solver is %s with a time of %.2f seconds.\n \n', solvers{idx}, min_time);
%% VERIFICATION
fprintf('VERIFICATION OF THE EQUILIBRIUM \n');

p.Lambda = (2*p.m_K*p.g)/(p.rho*p.S*p.v_w^2);
p.eps_c = p.c/p.l0;
[alpha_eq, ~, ~] = my_fzero("fun_equilibrio_paper",x_red0(3),1e-8,30,1e-6);

fprintf('Independent model results: \n theta_eq = %.4f rad = %.2f deg \n \n', alpha_eq, rad2deg(alpha_eq));


%% STABILITY WITH LINEAR THEORY
fprintf('STABILITY WITH LINEAR THEORY \n');

% Jacobian computation
dh = 1e-6;
J_eq = fun_jac_num(@RHS,0,x_eq,dh);
J_cond = cond(J_eq);
imag_J = any(imag(J_eq(:)) ~= 0);
if imag_J
    disp('Jacobian has imaginary values');
else
    disp('Jacobian has NO imaginary values');
end
fprintf('Jacobian condition number = %.2e \n', J_cond);

[vec_J, val_J]=eig(J_eq);

fprintf('Eigenvalues: \n')
eigenvalues_J = diag(val_J);

for i = 1:length(eigenvalues_J)
    n = real(eigenvalues_J(i));
    omega = imag(eigenvalues_J(i));
    if omega == 0
        fprintf('lambda_%d = %.2e\n', i, n);
    else
        fprintf('lambda_%d = %.2e + %.2ei\n', i, n, omega);
    end
end

%% INTEGRATOR
steps = 500;
tspan = linspace(1,100,steps);
opts = odeset('RelTol',1e-7,'AbsTol',1e-7);
[t,x] = ode15s(@RHS,tspan,x_eq,opts);

%% RESULTS
[F_S,M_OS,T_ASp,T_ASm,W_S,N,F_r,v_OS,ASp_AKp,ASm_AKm,up,um,ASp_OS,ASm_OS,OK_AKp,...
    OK_AKm,F_K,M_OK,v_OK,omega_KE,H_OK,W_K,T_AKp,T_AKm,F_a,M_a,alpha,beta,v_a,l_p,l_m,dxdt] = fun_get_results(x,p,frame);

%% PLOTS
figures = fun_draw_results_eq(x,t,p,frame);


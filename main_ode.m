%% MAIN FOR NUMERICAL INTEGRATION 
clc;
clear all;
close all;

% Define p as global
global p
global frame

% Import parameters
p = parameters;
frame = 0;

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


% v7 = fun_process_vector(vec_J(:,7))
% v8 = fun_process_vector(vec_J(:,8))
% v9 = fun_process_vector(vec_J(:,9))
% v10 = fun_process_vector(vec_J(:,10))
% v11 = fun_process_vector(vec_J(:,11))
% v12 = fun_process_vector(vec_J(:,12))
% v13 = fun_process_vector(vec_J(:,13))
% v14 = fun_process_vector(vec_J(:,14))
% v15 = fun_process_vector(vec_J(:,15))
% v16 = fun_process_vector(vec_J(:,16))
% v17 = fun_process_vector(vec_J(:,17))
% v18 = fun_process_vector(vec_J(:,18))



%% TRANSLATION OF OS
x_S_t = -4;
y_S_t = -6;
x_eq_t = zeros(18,1);

for i=1:18
    if i == 1
        x_eq_t(i) = x_S_t;
    elseif i == 2
        x_eq_t(i) = y_S_t;
    elseif i == 7
        x_eq_t(i) = x_eq(i) + x_S_t;
    elseif i == 8
        x_eq_t(i) = x_eq(i) + y_S_t;
    else
        x_eq_t(i) = x_eq(i);
    end
end
steps = 30;
tspan = linspace(1,60,steps);
opts = odeset('RelTol',1e-8,'AbsTol',1e-8);
[t,x] = ode15s(@RHS,tspan,x_eq_t,opts);

figures = fun_draw_results_eq(x,t,p,frame);
fun_download_fig('eq_translated',figures,360,350)
figures = 0;
%% ADDING A PERTURBATION
x_per = x_eq;

% Select which variable to perturb from 1 to 18 (bigger than than every
% component will be perturbed)
flag = 1;

for i = 1:18
    if flag <= 18 && flag == i
        if i == 3 || i == 10 || i == 11 || i == 12 || i == 16 || i == 17 || i == 18
            x_per(i) = x_eq(i) - 2*pi/360;
        else
            x_per(i) = x_eq(i) - 1e-3;
        end
    elseif flag > 18
        x_per = x_eq - 1e-3;
    end
end


% INTEGRATOR
steps = 300; 
tspan = linspace(1,300,steps);
opts = odeset('RelTol',1e-8,'AbsTol',1e-8);

[t,x] = ode15s(@RHS,tspan,x_per,opts);


% RESULTS
[F_S,M_OS,T_ASp,T_ASm,W_S,N,F_r,v_OS,ASp_AKp,ASm_AKm,up,um,ASp_OS,ASm_OS,OK_AKp,...
    OK_AKm,F_K,M_OK,v_OK,omega_KE,H_OK,W_K,T_AKp,T_AKm,F_a,M_a,alpha,beta,v_a,l_p,l_m,dxdt] = fun_get_results(x,p,frame);

% PLOTS
figures = fun_draw_results_stab(x, t, p, frame);
fun_download_fig('x_stab_xS',figures,360,350);
figures = 0;
%% NEUTRAL EIGENMODES
v15 = vec_J(:,15);
v17 = vec_J(:,17);
v18 = vec_J(:,18);

% Write the neutral eigenmode to be analyzed
v = v15;

x_n = x_eq + 1e-3*v;

steps = 700; 
tspan = linspace(1,100,steps);
opts = odeset('RelTol',1e-8,'AbsTol',1e-8);
[t,x] = ode15s(@RHS,tspan,x_n,opts);

figures = fun_draw_results_stab(x, t, p, frame);
fun_download_fig('x_neutral_15',figures,360,350);
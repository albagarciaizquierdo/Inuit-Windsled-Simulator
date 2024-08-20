clc;
clear all;
close all;

% Define p as global
global p
global frame

% Import parameters
p = parameters;
frame = 0;

%% EQUILIBRIUM STATE
% Initial guess for equilibrium
Gamma0 = 60*acos(-1)/180;
alpha0 = 10*acos(-1)/180;
x_K0 = - (p.l0+1)*cos(Gamma0);
z_K0 = - (p.l0+1)*sin(Gamma0);
x_red0 = [x_K0 z_K0 alpha0]';

% Equilibirum state calculation
[X_red_eq Error EXITO] = my_fzero("fun_equilibrio_red",x_red0,1e-8,30,1e-6);
x_eq = equilibrium_conditions(X_red_eq);

% Equilibrium state drawing
title_eq = "System in equilibrium";
box = 0;
fun_draw_system2(p,x_eq,title_eq,box)

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

format short e
eigenvalues_J = diag(val_J)

J_cond = cond(J_eq)
rho = zeros(1,18);
omega = zeros(1,18);

for i = 1:length(eigenvalues_J)
    rho(i) = real(eigenvalues_J(i));
    omega(i) = imag(eigenvalues_J(i));
end

% v7 = process_vector(vec_J(:,7))
% v8 = process_vector(vec_J(:,8))
% v9 = process_vector(vec_J(:,9))
% v10 = process_vector(vec_J(:,10))
% v11 = process_vector(vec_J(:,11))
% v12 = process_vector(vec_J(:,12))
% v13 = process_vector(vec_J(:,13))
% v14 = process_vector(vec_J(:,14))
% v15 = process_vector(vec_J(:,15))
% v16 = process_vector(vec_J(:,16))
% v17 = process_vector(vec_J(:,17))
% v18 = process_vector(vec_J(:,18))



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
% steps = 30;
% tspan = linspace(1,60,steps);
% opts = odeset('RelTol',1e-8,'AbsTol',1e-8);
% [t,x] = ode15s(@RHS,tspan,x_eq_t,opts);

% figures = fun_draw_results_eq(x,t,p,frame);
% fun_download_fig('eq_t',figures,360,350)

%% PERTURBATION
x_per = x_eq;

flag = 20;

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
%% INTEGRATOR

steps = 300; 
tspan = linspace(1,300,steps);
opts = odeset('RelTol',1e-8,'AbsTol',1e-8);
x0 = initial_conditions;
x0(4) = -1e-3;

[t,x] = ode15s(@RHS,tspan,x_per,opts);
% [t,x] = ode15s(@RHS,tspan,x0,opts);

%% RESULTS
[F_S,M_OS,T_ASp,T_ASm,W_S,N,F_r,v_OS,ASp_AKp,ASm_AKm,up,um,ASp_OS,ASm_OS,OK_AKp,...
    OK_AKm,F_K,M_OK,v_OK,omega_KE,H_OK,W_K,T_AKp,T_AKm,F_a,M_a,alpha,beta,v_a,l_p,l_m,dxdt] = fun_get_results(x,p,frame);

%% PLOTS
frame = 0; % SE (For SK frame = 1)
% figures = fun_draw_results(x,t,p,frame);
figures = fun_draw_results_stab(x, t, p, frame);
fun_download_fig('all_3',figures,360,400);

%% VIDEO
% fun_video(p,x,"xeq_50_fixed",t)

%% SAVE FIGURES
% % Create folder
% folderName = 'Figures_13_07_pert_vxS_SE_neg';
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
function output_vector = process_vector(input_vector)
    % Esta función recibe un vector de 18 filas, redondea a 0 los valores menores 
    % (en valor absoluto) de 1e-10, y convierte los demás números a notación 
    % científica con 3 cifras significativas, mostrando tanto la parte real como
    % la parte imaginaria.

    % Verificar que el input_vector tiene 18 filas
    if length(input_vector) ~= 18
        error('El vector de entrada debe tener 18 filas.');
    end
    
    % Inicializar el vector de salida
    output_vector = zeros(size(input_vector));
    
    % Procesar cada elemento del vector
    for i = 1:length(input_vector)
        real_part = real(input_vector(i));
        imag_part = imag(input_vector(i));
        
        % Redondear la parte real e imaginaria si son significativas
        if abs(real_part) < 1e-10
            real_part = 0;
        else
            real_part = round(real_part, 3, 'significant');
        end
        
        if abs(imag_part) < 1e-10
            imag_part = 0;
        else
            imag_part = round(imag_part, 3, 'significant');
        end
        
        % Reconstruir el número complejo con las partes redondeadas
        output_vector(i) = real_part + 1i * imag_part;
    end
    
end

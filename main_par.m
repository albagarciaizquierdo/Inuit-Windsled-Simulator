%% PARAMETRIC STUDY
clc;
clear all;
close all;

% Define p as global
global p
frame = 0; % SE

% Import parameters
p = parameters;

% Plot settings
color.blue = [0 0 1];
color.red = [1 0 0];
color.green = [0 0.6 0];
color.orange = [1 0.7 0.1];
color.purple = [0.5 0.2 0.6];

set(groot,'defaultAxesFontSize', 11);
set(groot,'defaultLineLineWidth', 1);
set(groot,'defaultAxesColorOrder', [color.blue; color.red; color.green; color.orange; color.purple]);
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
[X_red_eq, ~, ~] = my_fzero("fun_equilibrio_red",x_red0,1e-8,30,1e-6);
x_eq = equilibrium_conditions(X_red_eq);
[gamma_p_eq, gamma_m_eq] = fun_gamma(x_eq,p);

%% RESULTS
points = 100;
vw_max = 14;
fig = 0;

% Flag = 1 --> Kite position vs wind velocity with vw_min 
% figures = fun_par_fig(1,p,points,vw_max,frame,x_red0,color);

% Flag = 2 --> Length and angles for filtered velocity vs vw
% figures = fun_par_fig(2,p,points,vw_max,frame,x_red0,color);

% Flag = 3 --> All vs vw for different l0
% figures = fun_par_fig(3,p,points,vw_max,frame,x_red0,color);

% Flag = 4 --> z_K vs vw for different m_S
% figures = fun_par_fig(4,p,points,vw_max,frame,x_red0,color);

% Flag = 5 --> All vs vw for different S
figures = fun_par_fig(5,p,points,vw_max,frame,x_red0,color);

% Flag = 6 --> All vs vw for different yA_K
% figures = fun_par_fig(6,p,points,vw_max,frame,x_red0,color);

% Flag = 7 --> All vs vw for different zA_K
% figures = fun_par_fig(7,p,points,vw_max,frame,x_red0,color);

% Flag = 8 --> All vs vw for different xA_K
% figures = fun_par_fig(8,p,points,vw_max,frame,x_red0,color);

% Flag = 9 --> All vs vw for different w_S
% figures = fun_par_fig(9,p,points,vw_max,frame,x_red0,color);

% Flag = 10 --> All vs vw for different l_S
% figures = fun_par_fig(10,p,points,vw_max,frame,x_red0,color);

% Flag = 11 --> z_K vs vw for different xA_S
% figures = fun_par_fig(11,p,points,vw_max,frame,x_red0,color);

% Flag = 12 --> z_K vs vw for different yA_S
% figures = fun_par_fig(12,p,points,vw_max,frame,x_red0,color);

% Flag = 13 --> z_K vs vw for different mu_s
% figures = fun_par_fig(13,p,points,vw_max,frame,x_red0,color);

%% SAVE FIGURES
fig = fig + length(figures);

foldername = 'prueba';
width = 360;
height = 350;
fun_download_fig(foldername,figures,width,height)
